//! Isobaric labelling kits — TMT, TMTpro, iTRAQ, DiLeu — with every channel's reporter-ion m/z.
//!
//! | You want to | Call | mzLib type |
//! |---|---|---|
//! | List every kit mzLib can name, with its channels | [`kits`]`(None)` | `IsobaricMassTag`, `IsobaricMassTagType` |
//! | One kit's channels and reporter-ion m/z | [`kits`]`(Some("TMT18"))` | `IsobaricMassTag.TryGetTagType` |
//! | The channels grouped per kit | [`IsobaricKits::by_kit`], [`IsobaricKits::kit`] | `IsobaricMassTag.ChannelLabels`, `ReporterIonMzs` |
//! | The window a reporter intensity is read in | [`IsobaricKits::absolute_tolerance`], [`ReporterChannel::mz_min`] | `IsobaricMassTag.AbsoluteToleranceValue` |
//!
//! An isobaric experiment is read out of the low-m/z reporter ions: one ion per channel, a few
//! millidaltons apart. Getting a channel's label or its m/z wrong mislabels every sample after it,
//! so this crate keeps no table of its own. [`kits`] returns mzLib's:
//!
//! ```
//! # mzlib_replay::activate();
//! let tmtpro = mzlib::isobaric::kits(Some("TMT18"))?;
//! let kit = &tmtpro.by_kit()?[0];
//! assert_eq!((kit.name.as_str(), kit.channel_count), ("TMT18", 18));
//!
//! let labels: Vec<&str> = kit.channels[..3].iter().map(|c| c.label.as_str()).collect();
//! assert_eq!(labels, ["126", "127N", "127C"]);
//! assert!((kit.channels[0].reporter_ion_mz - 126.12773).abs() < 1e-5);    // m/z, charge 1
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **Nothing in that table is typed in.** mzLib holds only the channel labels. Each m/z is a
//! `DI HCD` diagnostic-ion line of the kit's `Multiplex Label` modification in mzLib's embedded
//! `TMT.txt`, plus one proton, sorted ascending and paired with the labels by position. The m/z are
//! therefore the same ones MetaMorpheus uses to quantify a TMT search. They are **theoretical**,
//! not calibrated or observed.
//!
//! The kit is looked up the way MetaMorpheus looks it up: by its whole name, case-insensitively,
//! optionally followed by `" on <motif>"` — `"TMT10"`, `"tmt10"`, `"TMT6-plex"`,
//! `"iTRAQ-4plex on K"`. Never by substring, so `"TMT10plex"` is refused rather than mistaken for a
//! kit whose name it contains.
//!
//! ## Every kit at once
//!
//! ```
//! # mzlib_replay::activate();
//! let catalogue = mzlib::isobaric::kits(None)?;
//! let names: Vec<(&str, u64)> = catalogue
//!     .kits
//!     .iter()
//!     .map(|k| (k.kit.as_str(), k.channel_count))
//!     .collect();
//! assert_eq!(names[..2], [("TMT6", 6), ("TMT10", 10)]);
//! assert_eq!(catalogue.absolute_tolerance, 0.003);                       // Da
//!
//! // iTRAQ 8-plex has no 120 channel: the phenylalanine immonium ion sits at 120.081.
//! assert_eq!(
//!     catalogue.kit("iTRAQ8")?.labels(),
//!     ["113", "114", "115", "116", "117", "118", "119", "121"]
//! );
//! # Ok::<(), mzlib::MzLibError>(())
//! ```
//!
//! **What the window means.** [`ReporterChannel::mz_min`] and [`ReporterChannel::mz_max`] are the
//! reporter m/z minus and plus [`IsobaricKits::absolute_tolerance`]: mzLib reads a reporter
//! intensity as the most intense peak inside that window. TMT16 is the lowest sixteen channels of
//! TMTpro 18-plex and has no modification entry of its own.
//!
//! **What this does not do.** It names the channels; it does not read reporter intensities out of
//! spectra, and it does not build an isobaric experimental design. An SDRF's label-free design is
//! [`crate::sdrf::design_with`].
//!
//! Ported from pyMzLib's `pymzlib.isobaric`, which decided the verb, the wire fields and the caveats.
//!
//! ## Cite
//!
#![doc = include_str!("../docs/reference/cite.isobaric.md")]

use std::time::Duration;

use serde::Deserialize;

use crate::bridge::{self, MzLibError, Result};
use crate::readers::Table;

/// The default timeout, matching pyMzLib's `isobaric.kits`. No file, no network.
const DEFAULT_TIMEOUT: Duration = Duration::from_secs(60);

/// One kit in [`IsobaricKits::kits`]: mzLib's name for it and how many channels it has.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct KitSummary {
    /// mzLib's `IsobaricMassTagType` member, e.g. `"TMT18"`. `TMT16` and `TMT18` are TMTpro.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub kit: String,
    /// How many channels the kit has.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub channel_count: u64,
}

/// One channel of a kit and the theoretical m/z of its reporter ion: a row of
/// [`IsobaricKits::channels`].
#[derive(Debug, Clone, PartialEq)]
pub struct ReporterChannel {
    /// mzLib's kit name (an `IsobaricMassTagType` member, e.g. `"TMT18"`).
    pub kit: String,
    /// 0-based position within the kit, ascending with `reporter_ion_mz`.
    pub index: u32,
    /// The channel's name as MetaMorpheus spells it — `"126"`, `"127N"`, `"115a"`.
    pub label: String,
    /// Theoretical m/z at charge 1: the diagnostic-ion mass from mzLib's `TMT.txt` plus one
    /// proton. Not a calibrated or observed value.
    pub reporter_ion_mz: f64,
    /// `reporter_ion_mz` minus the matching tolerance, in m/z.
    pub mz_min: f64,
    /// `reporter_ion_mz` plus the matching tolerance, in m/z.
    pub mz_max: f64,
}

/// One kit and its channels, in ascending reporter-ion m/z: an element of
/// [`IsobaricKits::by_kit`].
#[derive(Debug, Clone, PartialEq)]
pub struct IsobaricKit {
    /// mzLib's kit name, an `IsobaricMassTagType` member.
    pub name: String,
    /// How many channels the kit has.
    pub channel_count: u64,
    /// Every channel, in ascending reporter-ion m/z.
    pub channels: Vec<ReporterChannel>,
}

impl IsobaricKit {
    /// The channel labels, in ascending reporter-ion m/z.
    #[must_use]
    pub fn labels(&self) -> Vec<&str> {
        self.channels.iter().map(|c| c.label.as_str()).collect()
    }

    /// The reporter-ion m/z values, ascending, in m/z.
    #[must_use]
    pub fn reporter_ion_mzs(&self) -> Vec<f64> {
        self.channels.iter().map(|c| c.reporter_ion_mz).collect()
    }
}

/// The kits [`kits`] returned, as a long table and as structs.
///
/// The table has one row per (kit, channel), kits in mzLib's order and channels in ascending m/z:
/// `kit`, `channel_index`, `channel_label`, `reporter_ion_mz`, `mz_min` and `mz_max`, with the
/// meanings of [`ReporterChannel`] (the three m/z columns in m/z). [`Self::by_kit`] groups the
/// same rows by kit.
#[derive(Debug, Clone, Deserialize)]
pub struct IsobaricKits {
    /// The name you asked for, exactly as given; `None` when every kit was listed.
    #[serde(default)]
    pub kit: Option<String>,
    /// Kits listed.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub kit_count: u64,
    /// Channels listed, over every kit: the rows of [`Self::columns`].
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub record_count: u64,
    /// The half-width, in Da, of each channel's matching window
    /// (`IsobaricMassTag.AbsoluteToleranceValue`). mzLib reads a reporter intensity as the most
    /// intense peak within it.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub absolute_tolerance: f64,
    /// One summary per kit, in mzLib's `IsobaricMassTagType` order.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub kits: Vec<KitSummary>,
    /// The channel table.
    #[serde(flatten)]
    pub columns: Table,
    /// What the m/z are and are not — read these once.
    #[serde(default, deserialize_with = "bridge::null_to_default")]
    pub caveats: Vec<String>,
}

impl IsobaricKits {
    /// Every channel as a [`ReporterChannel`], kits in mzLib's order.
    ///
    /// # Errors
    ///
    /// [`MzLibError::Protocol`] if a column is missing or not the type the wire contract says.
    pub fn channels(&self) -> Result<Vec<ReporterChannel>> {
        if self.columns.names().is_empty() {
            return Ok(Vec::new());
        }
        let kit = self.columns.strings("kit")?;
        let index = self.columns.integers("channel_index")?;
        let label = self.columns.strings("channel_label")?;
        let mz = self.columns.floats("reporter_ion_mz")?;
        let mz_min = self.columns.floats("mz_min")?;
        let mz_max = self.columns.floats("mz_max")?;
        let missing = |what: &str, row: usize| {
            MzLibError::Protocol(format!("isobaric kits returned no {what} on row {row}"))
        };
        (0..kit.len())
            .map(|i| {
                Ok(ReporterChannel {
                    kit: kit[i].clone().ok_or_else(|| missing("kit", i))?,
                    index: index[i]
                        .and_then(|v| u32::try_from(v).ok())
                        .ok_or_else(|| missing("channel_index", i))?,
                    label: label[i]
                        .clone()
                        .ok_or_else(|| missing("channel_label", i))?,
                    reporter_ion_mz: mz[i].ok_or_else(|| missing("reporter_ion_mz", i))?,
                    mz_min: mz_min[i].ok_or_else(|| missing("mz_min", i))?,
                    mz_max: mz_max[i].ok_or_else(|| missing("mz_max", i))?,
                })
            })
            .collect()
    }

    /// Every kit as an [`IsobaricKit`], in mzLib's order.
    ///
    /// # Errors
    ///
    /// As [`Self::channels`].
    pub fn by_kit(&self) -> Result<Vec<IsobaricKit>> {
        let channels = self.channels()?;
        Ok(self
            .kits
            .iter()
            .map(|summary| IsobaricKit {
                name: summary.kit.clone(),
                channel_count: summary.channel_count,
                channels: channels
                    .iter()
                    .filter(|c| c.kit == summary.kit)
                    .cloned()
                    .collect(),
            })
            .collect())
    }

    /// The listed kit whose name is exactly `name` (mzLib's spelling, e.g. `"iTRAQ8"`).
    ///
    /// # Errors
    ///
    /// [`MzLibError::Usage`] if no listed kit has that name; otherwise as [`Self::channels`].
    pub fn kit(&self, name: &str) -> Result<IsobaricKit> {
        self.by_kit()?
            .into_iter()
            .find(|k| k.name == name)
            .ok_or_else(|| {
                let listed: Vec<&str> = self.kits.iter().map(|k| k.kit.as_str()).collect();
                MzLibError::Usage(format!(
                    "No kit named {name:?} in this result; it lists {listed:?}."
                ))
            })
    }
}

/// The isobaric kits mzLib can name, with every channel's label and reporter-ion m/z.
///
/// Calls mzLib's `IsobaricMassTag.TryGetIsobaricMassTag` for each `IsobaricMassTagType`, or for
/// the one `kit` resolves to through `IsobaricMassTag.TryGetTagType`. No file, no network.
///
/// `kit` is one kit, matched by mzLib's whole-name rule: case-insensitive, optionally followed by
/// `" on <motif>"` (`"TMT10"`, `"TMT6-plex"`, `"iTRAQ-4plex on K"`). `None` lists every kit.
#[doc = include_str!("../docs/reference/isobaric.kits.md")]
///
/// # Errors this crate adds
///
/// [`MzLibError::Usage`] for a blank `kit`, before anything is spawned.
///
/// # Examples
///
/// TMTpro 18-plex, whose top channel sits at 135.15:
///
/// ```
/// # mzlib_replay::activate();
/// let tmtpro = mzlib::isobaric::kits(Some("TMT18"))?;
/// assert_eq!((tmtpro.kit.as_deref(), tmtpro.record_count), (Some("TMT18"), 18));
///
/// let last = tmtpro.channels()?.pop().unwrap();
/// assert_eq!(last.label, "135N");
/// assert!((last.reporter_ion_mz - 135.1516).abs() < 1e-4);              // m/z
/// assert!((last.mz_max - last.mz_min - 0.006).abs() < 1e-9);            // ±0.003 Da
/// # Ok::<(), mzlib::MzLibError>(())
/// ```
#[doc = include_str!("../docs/reference/isobaric.kits.see-also.md")]
pub fn kits(kit: Option<&str>) -> Result<IsobaricKits> {
    let args = kit_args(kit)?;
    let data = bridge::invoke(&args, None, Some(DEFAULT_TIMEOUT))?;
    serde_json::from_value(data).map_err(|error| {
        MzLibError::Protocol(format!(
            "isobaric kits payload could not be interpreted: {error}"
        ))
    })
}

fn kit_args(kit: Option<&str>) -> Result<Vec<String>> {
    let mut args = vec!["isobaric".to_owned(), "kits".to_owned()];
    if let Some(kit) = kit {
        if kit.trim().is_empty() {
            return Err(MzLibError::Usage(
                "kit must be a kit name such as 'TMT10' or 'iTRAQ-4plex', or None for every kit."
                    .to_owned(),
            ));
        }
        args.push("--kit".to_owned());
        args.push(kit.to_owned());
    }
    Ok(args)
}

#[cfg(test)]
mod tests {
    //! Offline, against pyMzLib's recordings of the real bridge, shared byte for byte.

    use super::*;

    fn recorded(text: &str) -> IsobaricKits {
        let value: serde_json::Value = serde_json::from_str(text).unwrap();
        let data = match value.get("data") {
            Some(data) if value.get("ok").is_some() => data.clone(),
            _ => value,
        };
        serde_json::from_value(data).unwrap()
    }

    fn catalogue() -> IsobaricKits {
        recorded(include_str!("../tests/fixtures/isobaric_kits.json"))
    }

    #[test]
    fn every_kit_lists_its_channels_in_mzlib_order() {
        let all = catalogue();
        assert_eq!(all.kit, None);
        assert_eq!(all.kit_count, 9);
        let names: Vec<&str> = all.kits.iter().map(|k| k.kit.as_str()).collect();
        assert_eq!(
            names,
            ["TMT6", "TMT10", "TMT11", "TMT16", "TMT18", "iTRAQ4", "iTRAQ8", "diLeu4", "diLeu12"]
        );
        let total: u64 = all.kits.iter().map(|k| k.channel_count).sum();
        assert_eq!(total, all.record_count);
        assert_eq!(all.channels().unwrap().len() as u64, all.record_count);
    }

    #[test]
    fn by_kit_groups_the_table_without_losing_a_row() {
        let all = catalogue();
        for kit in all.by_kit().unwrap() {
            assert_eq!(kit.channels.len() as u64, kit.channel_count, "{}", kit.name);
            let mzs = kit.reporter_ion_mzs();
            assert!(mzs.windows(2).all(|w| w[0] < w[1]), "{} ascends", kit.name);
            let indices: Vec<u32> = kit.channels.iter().map(|c| c.index).collect();
            assert_eq!(indices, (0..kit.channels.len() as u32).collect::<Vec<_>>());
        }
    }

    #[test]
    fn the_window_is_the_wire_window_not_one_recomputed_here() {
        let tmtpro = recorded(include_str!("../tests/fixtures/isobaric_kits_TMT18.json"));
        let first = &tmtpro.channels().unwrap()[0];
        assert_eq!(first.label, "126");
        // The wire's own values (serde_json's default float parse may differ in the last ulp).
        assert!((first.mz_min - 126.124_725_954_759_01).abs() < 1e-12);
        assert!((first.mz_max - 126.130_725_954_759_01).abs() < 1e-12);
        assert_eq!(tmtpro.absolute_tolerance, 0.003);
    }

    #[test]
    fn an_unlisted_kit_is_a_usage_error_naming_the_kits() {
        let error = catalogue().kit("TMT10plex").unwrap_err();
        assert!(
            matches!(error, MzLibError::Usage(ref m) if m.contains("TMT10")),
            "{error}"
        );
    }

    #[test]
    fn a_kit_reaches_the_bridge_and_a_blank_one_does_not() {
        assert_eq!(kit_args(None).unwrap(), ["isobaric", "kits"]);
        assert_eq!(
            kit_args(Some("iTRAQ-4plex on K")).unwrap(),
            ["isobaric", "kits", "--kit", "iTRAQ-4plex on K"]
        );
        assert!(matches!(kit_args(Some("  ")), Err(MzLibError::Usage(_))));
    }
}
