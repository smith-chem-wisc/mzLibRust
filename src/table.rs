//! A columnar table with typed accessors, shared by every verb that returns one.
//!
//! The column set of a result is not knowable at compile time — it depends on the file format for
//! [`crate::readers`] and on the model family for [`crate::prediction`] — so those verbs cannot
//! return a struct with named fields. They return this, whose accessors do the projection
//! properly: every cell is an [`Option`], a wire `null` becomes [`None`], and a column whose
//! values are not the type you asked for is an error rather than a silent default.

use std::collections::BTreeMap;

use serde_json::Value;

use crate::bridge::{MzLibError, Result};

/// A columnar table: one array per field, with typed accessors.
///
/// The column set depends on the format — a TopPIC file has 36 columns and a Crux file 23 — so
/// this cannot be a struct with named fields. What it can be is a map whose accessors do the
/// projection properly: every cell is an [`Option`], a wire `null` becomes [`None`], and a column
/// whose values are not the type you asked for is an error rather than a silent default.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Table {
    names: Vec<String>,
    columns: BTreeMap<String, Vec<Value>>,
}

impl Table {
    /// The column names, **in the order mzLib declares them** — base-class fields first, then
    /// each format's own.
    #[must_use]
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// The number of rows, or zero when the table went to disk instead.
    #[must_use]
    pub fn rows(&self) -> usize {
        self.names
            .first()
            .and_then(|name| self.columns.get(name))
            .map_or(0, Vec::len)
    }

    /// Whether the table carries any rows.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.rows() == 0
    }

    /// Whether a column is present.
    #[must_use]
    pub fn has(&self, name: &str) -> bool {
        self.columns.contains_key(name)
    }

    /// A column's raw JSON values, for a shape the typed accessors do not cover.
    #[must_use]
    pub fn raw(&self, name: &str) -> Option<&[Value]> {
        self.columns.get(name).map(Vec::as_slice)
    }

    /// A column as floating-point numbers, with `null` as [`None`].
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if the column is absent, [`MzLibError::Protocol`] if a value is
    /// present but is not a number.
    pub fn floats(&self, name: &str) -> Result<Vec<Option<f64>>> {
        self.project(name, "a number", |value| match value {
            Value::Null => Some(None),
            Value::Number(number) => number.as_f64().map(Some),
            _ => None,
        })
    }

    /// A column as whole numbers, with `null` as [`None`].
    ///
    /// # Errors
    ///
    /// As [`Table::floats`], and also when a value is a number with a fractional part — silently
    /// truncating one would turn a protocol change into wrong data.
    pub fn integers(&self, name: &str) -> Result<Vec<Option<i64>>> {
        self.project(name, "a whole number", |value| match value {
            Value::Null => Some(None),
            Value::Number(number) => number.as_i64().map(Some),
            _ => None,
        })
    }

    /// A column as strings, with `null` as [`None`].
    ///
    /// # Errors
    ///
    /// As [`Table::floats`], for a value that is not a string.
    pub fn strings(&self, name: &str) -> Result<Vec<Option<String>>> {
        self.project(name, "a string", |value| match value {
            Value::Null => Some(None),
            Value::String(text) => Some(Some(text.clone())),
            _ => None,
        })
    }

    /// A column as booleans, with `null` as [`None`].
    ///
    /// `None` genuinely means *unknown* in this library rather than *false* — Casanovo's
    /// `is_decoy` is the case that matters, because de novo sequencing has no target/decoy label
    /// at all and a `false` there would be a fabricated value someone could filter on.
    ///
    /// # Errors
    ///
    /// As [`Table::floats`], for a value that is not a boolean.
    pub fn booleans(&self, name: &str) -> Result<Vec<Option<bool>>> {
        self.project(name, "a boolean", |value| match value {
            Value::Null => Some(None),
            Value::Bool(flag) => Some(Some(*flag)),
            _ => None,
        })
    }

    /// A column whose every cell is itself an array of numbers.
    ///
    /// The shape [`crate::readers::read_spectra`] returns for `mz` and `intensity` under
    /// [`crate::readers::SpectraOptions::peaks`]: one array per scan, not one number per scan.
    ///
    /// # Errors
    ///
    /// As [`Table::floats`], for a value that is not an array of numbers.
    pub fn float_arrays(&self, name: &str) -> Result<Vec<Option<Vec<Option<f64>>>>> {
        self.project(name, "an array of numbers", |value| match value {
            Value::Null => Some(None),
            Value::Array(items) => items
                .iter()
                .map(|item| match item {
                    Value::Null => Some(None),
                    Value::Number(number) => number.as_f64().map(Some),
                    _ => None,
                })
                .collect::<Option<Vec<_>>>()
                .map(Some),
            _ => None,
        })
    }

    /// A column whose every cell is itself an array of strings.
    ///
    /// The shape a fragment prediction returns for `fragment_annotations`: one array per peptide,
    /// index-aligned with that row's `fragment_mz` and `fragment_intensity`.
    ///
    /// # Errors
    ///
    /// As [`Table::floats`], for a value that is not an array of strings.
    pub fn string_arrays(&self, name: &str) -> Result<Vec<Option<Vec<Option<String>>>>> {
        self.project(name, "an array of strings", |value| match value {
            Value::Null => Some(None),
            Value::Array(items) => items
                .iter()
                .map(|item| match item {
                    Value::Null => Some(None),
                    Value::String(text) => Some(Some(text.clone())),
                    _ => None,
                })
                .collect::<Option<Vec<_>>>()
                .map(Some),
            _ => None,
        })
    }

    /// The shared body of every typed accessor: absent is a usage error, wrong-typed is a protocol
    /// error, and the two are never conflated.
    fn project<T>(
        &self,
        name: &str,
        expected: &str,
        convert: impl Fn(&Value) -> Option<T>,
    ) -> Result<Vec<T>> {
        let column = self.columns.get(name).ok_or_else(|| {
            // Naming what IS there, because the column set is per-format and a caller who guessed
            // a name from another format has no other way to find out.
            MzLibError::Usage(format!(
                "No column '{name}' in this table. Its columns are: {}.",
                self.names.join(", ")
            ))
        })?;

        column
            .iter()
            .enumerate()
            .map(|(row, value)| {
                convert(value).ok_or_else(|| {
                    MzLibError::Protocol(format!(
                        "Column '{name}' row {row} is not {expected}: {value}"
                    ))
                })
            })
            .collect()
    }

    /// Builds a table from the wire's `column_names` and `columns`.
    ///
    /// Order comes from `column_names` rather than from the map, because a JSON object has no
    /// order a parser must preserve — reading it back from the map would give alphabetical
    /// columns, which is not the order mzLib declares its fields in.
    pub(crate) fn from_wire(
        names: Vec<String>,
        columns: Option<BTreeMap<String, Vec<Value>>>,
    ) -> Self {
        let columns = columns.unwrap_or_default();
        // Only names that actually have a column, so `names()` cannot promise a column `raw()`
        // then fails to return.
        let names = if names.is_empty() {
            columns.keys().cloned().collect()
        } else {
            names
                .into_iter()
                .filter(|name| columns.contains_key(name))
                .collect()
        };
        Self { names, columns }
    }
}
