# cherentrace

Modifications for CORSIKA, the IACT interface, and sim_telarray to enable
tracing back Cherenkov photons that reached the camera to their emission
points, and to more easily correlate them with particles that reach the ground.
Designed for use with ctapipe.

Requirements:
- CORSIKA 7.7500.
- IACT/ATMO (bernlohr) 1.67.
- sim_telarray 2024-02-28.
- pyeventio at least v1.15.0 (2025-03-06).
- ctapipe at least 0.23.1.

## Changes made by this package

- Instead of writing MASS and CHARGE in the emitter table, write PARTICLE_ID and
  GENERATION. These can be used for matches to the particle table.
- Instead of writing X and Y in the emitter table, write XEM and YEM. This is
  the emission position in the shower frame.
- Write an additional MC_PHOTONS block out from sim_telarray that records the
  Cherenkov photons that survived raytracing to the focal plane (plane of the
  camera). This block includes the emitter information of each photon if
  available.

A reminder of how exactly data gets from CORSIKA to a simtel file:
```
CORSIKA -> IACT interface -> eventio-format binary stream -> multipipe_corsika ---> write to file
                                                                                \-> sim_telarray with telescope class A -> eventio-format to file
                                                                                \-> sim_telarray with telescope class B -> eventio-format to file
                                                                                \-> ...
```
The "IACT interface" is a set of C routines called from Fortran, which convert
data into eventio format and pass it on to (generally) `multipipe_corsika`.

`multipipe_corsika` takes the eventio-format binary stream from the IACT
interface and splits it between multiple output locations, including
potentially saving directly to disk and being passed to sim_telarray to
simulate different types of telescopes or NSB levels from the same air shower.

sim_telarray writes data in eventio format to simtel files, but it doesn't
necessarily pass through all of the eventio data that it originally received:
some of it is reformatted into different data structures, and some is dropped.

## Usage

There are two "stages" of configuration of CORSIKA/IACT relevant to these tools:
CORSIKA compile flags and CORSIKA steering commands (keywords) issued in the
input card.

CORSIKA compile flags can be managed by the `build_all` script distributed with
sim_telarray. This is recommended as the easiest way (see next section).

Useful references:
- Chapters 3-4 in the [CORSIKA manual](https://web.iap.kit.edu/corsika/usersguide/usersguide.pdf).
- Chapter 3 in the [sim_telarray manual](https://www.mpi-hd.mpg.de/hfm/~bernlohr/sim_telarray/Documentation/sim_telarray.pdf).

Follow these steps:
1. **Patch:**
    ```sh
    patch --backup /path/to/install/dir/corsika-77500/src/corsika.F patches/corsika.F.patch
    patch --backup /path/to/install/dir/corsika-77500/bernlohr/iact.c patches/iact.c.patch
    patch --backup /path/to/install/dir/sim_telarray/common/sim_telarray.c patches/sim_telarray.c.patch
    ```

2. **Compile**:
    Run
    ```sh
    TAR_KEEP_NEWER=1 ./build_all prod6 qgs2 iactext muprod store-emitter
    ```

    The effect of each of these compilation flags is as follows:
    - `IACTEXT`: Pass photon "emitter" information from CORSIKA into the IACT
      interface.
    - `MUPROD`: Include "fated" muons (that don't survive to the observation
      level) in the particle table passed to the IACT interface.
    - `STORE_EMITTER`: Set the default value of the `IACT STORE-EMITTER` keyword
      to true.

    Fated muons are called "decaying" muons in official documentation. We use
    this more general term because the muons could have been removed from the
    particle stack for other reasons besides decay (recorded by the fate
    index).

3. **Simulate:**
    Add your selection of the following commands to your CORSIKA input card:
    - `MUADDI YES`: Include "birth" muon entries (point at which a muon was
      created) in the particle table passed to the IACT interface.
    - `IACT STORE-EMITTER YES`: Pass photon "emitter" information from the IACT
      interface into the eventio output stream. Defaults to `YES` if
      sim_telarray was compiled with the `STORE_EMITTER` flag.
    - `IACT STORE-PARTICLES YES`: Pass the particle table from the IACT
      interface into the eventio output stream.

    Another relevant command is `OBSLEV` which defines a layer
    (observation level) at which to record particle information. There can be
    up to 10 observation levels, listed in any order. The lowest observation
    level is taken to be the ground level that the telescopes are on. The
    counting of the level numbers (only relevant to the particle table) is such
    that level 1 is the highest level and level N is the lowest level.

    You can use the demo notebook to test your created simtel file.

4. **Analyse:**
    The Python module `cherentrace` provides two functions for accessing the
    photon and particle tables in a convenient format. These functions take care
    of transforming from CORSIKA units (cm) to ctapipe default units (m) and
    deriving other useful quantities.

    -
        ```python
        import cherentrace
        df_photons = cherentrace.get_photons(source, event, tel_id, to_telescope_frame = True)
        ```
        Get a Pandas DataFrame of the Cherenkov photons that reached the focal
        plane, including emitter information of each photon if available.

    -
        ```python
        import cherentrace
        df_particles = cherentrace.get_particles(source, event)
        ```
        Get a Pandas DataFrame of the particles that reached each observation
        level, including additional information about muons if available.


## Columns of the photon and particle tables

Note that all positions, vectors, and angles are reported in the CORSIKA
coordinate system: x is positive towards magnetic north, y is positive towards
west, and z is positive upwards. Azimuth is measured anti-clockwise from the x
axis (north).

### Photon table

|     Column    | Description |
|---------------|-------------|
| x             | Position in the camera frame (m or deg, depending on value of `to_telescope_frame`). |
| y             | See above. |
| alt           | Reconstructed arrival direction (deg) [only if `to_telescope_frame=True`]. This is calculated from the position on the camera using the focal length and may be affected by the PSF. |
| az            | See above. |
| cx            | Unit momentum vector at the camera (points back towards the mirror). |
| cy            | See above. |
| xem           | Emission x position with respect to shower core at ground level (m). |
| yem           | Emission y position with respect to shower core at ground level (m). |
| zem           | Emission z position with respect to sea level (m). |
| time          | Arrival time at ground relative to time when the primary travelling at v = c would arrive at the core in the CORSIKA detection plane (ns). |
| pixel_id      | Pixel ID if the photon was registered by a pixel, or -1 otherwise. |
| wavelength    | Cherenkov wavelength. |
| particle_id   | Emitting particle ID (see [Table 4 in the CORSIKA manual](https://web.iap.kit.edu/corsika/usersguide/usersguide.pdf#page=132) and [Particle ID](#particle-id) below). |
| generation    | Emitting particle's complete generation counter, see [Generation counter](#generation-counter) below. |
| emission_time | Emission time, counting since the primary entered the atmosphere or the first interaction (ns). |
| energy        | Emitting particle energy (GeV). |
| is_muon       | Convenience boolean that is True if the emitting particle was a muon, False otherwise. |


### Particle table

|    Column   | Description |
|-------------|-------------|
| x           | x position with respect to shower core at ground level (m). |
| y           | y position with respect to shower core at ground level (m). |
| z           | z position with respect to sea level (m). |
| cx          | Unit momentum vector at the reported position. |
| cy          | See above. |
| cz          | See above. |
| time        | Time since the primary entered the atmosphere or the first interaction (ns). |
| momentum    | Momentum (GeV/c). |
| weight      | The thinning weight or 1.0. |
| particle_id | Particle ID (see [Table 4 in the CORSIKA manual](https://web.iap.kit.edu/corsika/usersguide/usersguide.pdf#page=132) and [Particle ID](#particle-id) below). |
| obs_level   | Numbered observation level that the particle track passed through or -1 if this is an additional muon information entry. Note that the highest-numbered observation level is the ground level. |
| fate_index  | Reason that a muon track ended or -1, see [Fate index](#fate-index) below. |
| generation  | Last two or three digits of the generation counter, see [Generation counter](#generation-counter) below. |
| is_muon     | Convenience boolean that is True if this is any type of muon entry, False otherwise. |

### Particle ID

References:
- Table 4 in the [CORSIKA manual](https://web.iap.kit.edu/corsika/usersguide/usersguide.pdf#page=132)
- The [MUPROD option manual](https://web.iap.kit.edu/heck/publications/kit-swp-5_muprod.pdf)

The emitter information in the photon table only identifies muons as ID 5 or 6
(μ+ and μ-, respectively). In the particle table, there may be additional muon
information with other IDs.

The following table lists the meaning of each muon-related particle ID in the
particle table.

|   ID  | Meaning                                                                        |
|-------|--------------------------------------------------------------------------------|
|  5/6  | Point where a surviving muon track passed through the observation level.       |
| 75/76 | Point where a surviving muon track started (muon birth point).                 |
| 85/86 | Point where a fated muon track started (muon birth point).                     |
| 95/96 | Point where a fated muon track ended (decay, fatal interaction, or E/ang cut). |

The following table lists which types of muon entries are written to the output
for each combination of compilation and steering options.

| Steering command 🠋 / Compile flag 🠊 | None        | `MUPROD`                    |
|--------------------------------------|-------------|-----------------------------|
| None                                 | 5/6         | 5/6 + 95/96                 |
| `MUADDI YES`                         | 5/6 + 75/76 | 5/6 + 75/76 + 85/86 + 95/96 |

In other words, `MUPROD` enables information about muons that don't survive to
the next observation level. `MUADDI YES` enables information about muon birth
points.

### Fate index

Reference:
- Section 3.1 in the [MUPROD option manual](https://web.iap.kit.edu/heck/publications/kit-swp-5_muprod.pdf#page=11)

| Index | Meaning                                               |
|-------|-------------------------------------------------------|
|   1   | Muon track ends because of decay.                     |
|   2   | Muon track ends because of fatal nuclear interaction. |
|   3   | Muon track ends because of energy or angular cut.     |

### Generation counter

The generation counter tracks the interaction history of each particle by adding
specific units together to represent different types of interactions. In the
photon table we get the full generation counter, however in the particle table
the generation counter is only the last two or three digits of the generation
counter (two digits if an observation level is given, three digits if the
observation level is -1).

For a reference on the full generation counter, refer to Section 3.5.15 in the
[CORSIKA manual](https://web.iap.kit.edu/corsika/usersguide/usersguide.pdf#page=63).
Although this section is describing the `EHISTORY` option, this is still the
relevant description of the full generation counter that we get in the photon
table.

For a reference on how hadronic interactions are handled in the generation
counter, refer to Section 4 in the
[MUPROD option manual](https://web.iap.kit.edu/heck/publications/kit-swp-5_muprod.pdf#page=13).

Note that the last three digits of the generation counter are for the hadronic
interactions, so in the particle table (which only lists the last two or three
digits), only hadronic interactions are seen.
