# cherentrace

Modifications for CORSIKA, the IACT interface, and sim_telarray to enable
tracing back Cherenkov photons that reached the camera to their emission
points, and to more easily correlate them with particles that reach the ground.
Designed for use with ctapipe.

## Changes made by this package

- Instead of writing MASS and CHARGE in the emitter table, write PARTICLE_ID and
  GENERATION. These can used for matches to the particle table.
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
`multipipe_corsika` takes the eventio-format binary stream from the IACT
interface and splits it between multiple output locations, including
potentially saving directly to disk and being passed to sim_telarray to
simulate different types of telescopes or NSB levels from the same air shower.

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
    todo

2. **Compile**:
    Run
    ```sh
    TAR_KEEP_NEWER=1 ./build_all prod6 qgs2 iactext muprod store-emitter
    ```

    The effect of each of these compilation flags is as follows:
    - `IACTEXT`: Pass photon "emitter" information from CORSIKA into the IACT
      interface.
    - `MUPROD`: Include "fated" muons (that don't survive to the observation level)
      in the particle table passed to the IACT interface.
    - `STORE_EMITTER`: Set the default value of the `IACT STORE-EMITTER` keyword to
      true.

    Fated muons are called "decaying" muons in official documentation. We use this
    more general term because the muons could have been removed from the particle
    stack for other reasons besides decay (recorded by the fate index).

3. **Simulate:**
    Add your selection of the following commands to your CORSIKA input card:
    - `MUADDI YES`: Include "birth" muon entries (point at which a muon was created)
      in the particle table passed to the IACT interface.
    - `IACT STORE-EMITTER YES`: Pass photon "emitter" information from the IACT
      interface into the eventio output stream. Defaults to `YES` if sim_telarray
      was compiled with the `STORE_EMITTER` flag.
    - `IACT STORE-PARTICLES YES`: Pass the particle table from the IACT interface
      into the eventio output stream.

    Another relevant command is `OBSLEV` which defines a layer (observation level)
    at which to record particle information. There can be up to 10 observation
    levels, listed in any order. The lowest observation level is taken to be the
    ground level that the telescopes are on. The counting of the level numbers
    (only relevant to the particle table) is such that level 1 is the highest level
    and level N is the lowest level.

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
        Get a Pandas DataFrame of the Cherenkov photons that reached the focal plane,
        including emitter information of each photon if available.

    -
        ```python
        import cherentrace
        df_particles = cherentrace.get_particles(source, event)
        ```
        Get a Pandas DataFrame of the particles that reached each observation level,
        including additional information about muons if available.


## Columns of the photon and particle tables

### A note about Particle IDs

References:
- Table 4 in the [CORSIKA manual](https://web.iap.kit.edu/corsika/usersguide/usersguide.pdf#page=132)
- The [MUPROD option manual](https://web.iap.kit.edu/heck/publications/kit-swp-5_muprod.pdf)

The emitter information in the photon table only identifies muons as ID 5 or 6
(μ+ and μ-, respectively). In the particle table, there may be additional muon
information with other IDs.

The following table lists the meaning of each particle ID in the particle table.

|   ID  | Meaning                                                                        |
|-------|--------------------------------------------------------------------------------|
|  5/6  | Point where a surviving muon track passed through the observation level.       |
| 75/76 | Point where a surviving muon track started (muon birth point).                 |
| 85/86 | Point where a fated muon track started (muon birth point).                     |
| 95/96 | Point where a fated muon track ended (decay, fatal interaction, or E/Θ cut).   |

The following table lists which types of muon entries are written to the output
for each combination of compilation and steering options.

| Steering command 🠋 / Compile flag 🠊 | None        | `MUPROD`                    |
|--------------------------------------|-------------|-----------------------------|
| None                                 | 5/6         | 5/6 + 95/96                 |
| `MUADDI YES`                         | 5/6 + 75/76 | 5/6 + 75/76 + 85/86 + 95/96 |

In other words, `MUPROD` enables information about muons that don't survive to
the next observation level. `MUADDI YES` enables information about muon birth
points.

### Photon table

| Column | Value |
|--------|-------|
| todo   | todo  |

### Particle table

| Column | Value |
|--------|-------|
| todo   | todo  |
