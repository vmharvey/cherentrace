import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.time import Time
from ctapipe.coordinates import TelescopeFrame


def _assert_event_match(source, event):
  if (source.file_.header['run'] != event.index.obs_id
    or source.file_.current_mc_event_id != event.index.event_id):
    raise RuntimeError("EventSource is currently set to a different event than "
      "the provided event. This function should be called from within the loop "
      "over events from EventSource")
  else:
    pass

def _find_paired_data(df, particle_id):
  """Birth muon rows in the particle table (75, 76, 85, and 86) have a matching
  row immediately following them. This function looks up those rows and returns
  them. Sometimes, a few rows won't have a matching pair (possible CORSIKA
  bug?); this function does the extra checking necessary to confirm the match
  is correct. Indices of rows that don't have a match are also returned.

  Returns:
    keep_idx :  Pandas Index of the birth muons that have a match
    keep_vals:  Pandas DataFrame of the matching muon data for the birth muons
    drop_idx :  Pandas Index of the birth muons that do not have a match
  """
  if particle_id in [75, 76]:
    expected_pid = particle_id - 70
  elif particle_id in [85, 86]:
    expected_pid = particle_id + 10
  else:
    raise RuntimeError(f"Unexpected particle ID: {particle_id}")

  idx = df.query(f"particle_id == {particle_id}").index
  paired_idx = idx + 1
  found_idx = paired_idx.intersection(df.index)

  found_pid = df.loc[found_idx, 'particle_id']
  match_mask = found_pid == expected_pid

  keep_idx = found_idx[match_mask] - 1
  keep_vals = df.loc[found_idx[match_mask], :]
  drop_idx = found_idx[~match_mask]

  # Due to unknown reasons (maybe a CORSIKA bug?) sometimes the expected
  # matching entry in the next row is not found in the table, or it has the
  # wrong particle ID. There are usually only a handful of these and it
  # shouldn't be a big problem
  miss_idx = paired_idx.difference(df.index)
  drop_idx = drop_idx.union(miss_idx)
  if not drop_idx.empty:
    print(
      "Warning: No matching row found for muon entries with "
      f"ID={particle_id}: {drop_idx.values}"
    )

  return keep_idx, keep_vals, drop_idx

def _get_obslev_offsets(source, event):
  obslevs = get_corsika_obslevs(source)
  simshower = event.simulation.shower

  tanth = np.tan(90*u.deg - simshower.alt)
  xoff = -1*(obslevs - obslevs[-1])*tanth*np.cos(simshower.az) # not sure why this needs the minus sign
  yoff = (obslevs - obslevs[-1])*tanth*np.sin(simshower.az)

  # Cast to float32 to match the precision of the particle table x/y coordinates
  return {'x': xoff.value.astype(np.float32), 'y': yoff.value.astype(np.float32)}

def get_corsika_obslevs(source):
  input_card = source.file_.corsika_inputcards[0].decode('utf8')
  input_card_wordlines = [[w for w in l.split(' ') if w] for l in input_card.split('\n') if l and l[0] != '*']
  obslevs = sorted([float(l[1])/100 for l in input_card_wordlines if l and l[0] == 'OBSLEV'])
  # OBSLEV 1 must be the highest level, OBSLEV n must be the lowest level
  # Cast to float32 to match the precision of the particle table z coordinate
  return np.flip(np.asarray(obslevs, dtype = np.float32))

def get_photons(source, event, tel_id, to_telescope_frame = True):
  """Return a Pandas DataFrame containing the Cherenkov photons that reached the
  camera plane. If CORSIKA was compiled with the IACTEXT option and
  sim_telarray was compiled with the store-emitter option (todo check this), then
  emitter information for each photon will be included in the table.

  If to_telescope_frame is True, the photon x/y coordinates will be
  transformed into the telescope frame, and their absolute alt/az arrival
  directions will also be inserted into the DataFrame. If to_telescope_frame is
  False, the coordinates are left in the camera frame.

  Note this function must be used within the loop over ArrayEventContainers
  from the EventSource, to find the correct data in the underlying simtel file.

  Photon emission point coordinates are in the shower frame. Add the event's
  core position to offset the photon coordinates into the array frame.
  """

  _assert_event_match(source, event)

  try:
    true_photons = source.file_.current_photons[tel_id - 1]
  except KeyError:
    raise RuntimeError("No photon data found")

  true_emitter = None
  try:
    true_emitter = source.file_.current_emitter[tel_id - 1]
  except:
    pass

  df_photons = pd.DataFrame(true_photons)
  df_photons.rename(columns = {'photons': 'pixel_id'}, inplace = True)

  df_photons.x = df_photons.x/100 # convert to m
  df_photons.y = df_photons.y/100
  df_photons.zem = df_photons.zem/100
  df_photons.pixel_id = df_photons.pixel_id.astype(int)

  if to_telescope_frame:
    cam_frame = source.subarray.tels[tel_id].camera.geometry.frame
    tel_frame = TelescopeFrame(telescope_pointing = SkyCoord(
      alt = event.pointing.tel[tel_id].altitude,
      az = event.pointing.tel[tel_id].azimuth,
      frame = AltAz(
        obstime = Time.now(),
        location = EarthLocation.of_site('Roque de los Muchachos'),
      ),
    ))

    photon_x = u.Quantity(df_photons.x, u.m)
    photon_y = u.Quantity(df_photons.y, u.m)
    coord = SkyCoord(photon_x, photon_y, frame = cam_frame)
    trans = coord.transform_to(tel_frame)
    df_photons.x = trans.fov_lon.to_value('deg')
    df_photons.y = trans.fov_lat.to_value('deg')
    df_photons['alt'] = event.pointing.tel[tel_id].altitude.to_value('deg') + df_photons.y.values
    df_photons['az'] = event.pointing.tel[tel_id].azimuth.to_value('deg') + df_photons.x.values

  if true_emitter is not None:
    df_emitter = pd.DataFrame(true_emitter)
    # Unused
    df_emitter.drop(columns = ['time', 'wavelength'], inplace = True)
    # Edits to CORSIKA/IACT and sim_telarray mean these will be emission points
    df_emitter.rename(columns = {'x': 'xem', 'y': 'yem'}, inplace = True)
    df_emitter.xem = df_emitter.xem/100 # convert to m
    df_emitter.yem = df_emitter.yem/100
    df_emitter.emission_time = df_emitter.emission_time*1e9 # convert to ns
    # Edits to CORSIKA/IACT and sim_telarray mean these will be ID and generation
    df_emitter.rename(columns = {'mass': 'particle_id', 'charge': 'generation'}, inplace = True)
    df_emitter.particle_id = df_emitter.particle_id.astype(int)
    df_emitter.generation = df_emitter.generation.astype(int)

    df_emitter['is_muon'] = False
    df_emitter.loc[df_emitter.query("particle_id in [5, 6]").index, 'is_muon'] = True

    # Zip emitter data with photon data
    df_photons = pd.concat([df_photons, df_emitter], axis = 'columns')

  if true_emitter is None:
    new_col_order = [0, 1, 8, 9, 2, 3, 5, 4, 6, 7]
  else:
    new_col_order = [0, 1, 8, 9, 2, 3, 10, 11, 5, 4, 6, 7, 12, 13, 14, 15, 16]

  df_photons = df_photons[ df_photons.columns[new_col_order] ]

  return df_photons

def get_particles(source, event):
  """Return a Pandas DataFrame containing the particles that passed through
  each CORSIKA observation level. CORSIKA must be run with "IACT STORE-PARTICLES YES".
  If CORSIKA was run with "MUADDI YES" and/or compiled with the "MUPROD" option,
  the table will contain additional information about muons stored in particles
  of ID 75-76, 85-86, or 95-96, depending on the type of additional information
  that it is.

  Particle coordinates are in the shower frame. Add the event's core position to
  offset the particle coordinates into the array frame.
  """

  _assert_event_match(source, event)

  try:
    obslev_particles = source.file_.current_obslev_particles
  except AttributeError:
    raise RuntimeError("eventio version is too old, could not read the "
      "obslev_particles block")

  # CORSIKA must be run with "IACT STORE-PARTICLES YES"
  if obslev_particles is None:
    return None

  obslev_z = get_corsika_obslevs(source)
  obslev_xy = _get_obslev_offsets(source, event)

  df = pd.DataFrame(obslev_particles)
  df.drop_duplicates( # Remove duplicate 75/76 entries, which can also overlap with 85/86
    subset = ['cx', 'cy', 'momentum', 'time'],
    keep = 'first', inplace = True)
  df.x = df.x/100 # convert to m
  df.y = df.y/100
  df['cz'] = -1*np.sqrt(1 - df.cx**2 - df.cy**2) # downwards

  df['z'] = np.nan
  df['obs_level'] = -1
  df['fate_index'] = -1
  df['generation'] = (df.particle_id % 1000).astype(int)
  df.particle_id = ((df.particle_id - df.generation) / 1000).astype(int)

  # TODO: Add support for ID<0, which will be EHISTORY mother and grandmother
  # particles

  # For 0<ID<75,ID>100: gen number g, obs level number l: g×10 + l
  std_part = df.query("(particle_id > 0 and particle_id < 75) or particle_id > 100")
  _obs_level = std_part.generation % 10
  _generation = (std_part.generation - _obs_level) / 10 # Must be before the next line
  _obs_level = _obs_level.mask(_obs_level == 0, other = 10) # SPECIAL CASE OF 10 LEVELS
  df.loc[std_part.index, 'obs_level'] = _obs_level.astype(int)
  df.loc[std_part.index, 'generation'] = _generation.astype(int)
  df.loc[std_part.index, 'x'] = std_part.x - obslev_xy['x'][_obs_level - 1]
  df.loc[std_part.index, 'y'] = std_part.y - obslev_xy['y'][_obs_level - 1]
  df.loc[std_part.index, 'z'] = obslev_z[_obs_level - 1]

  addi_muon = df.query("75 <= particle_id <= 96")
  df.loc[addi_muon.index, 'z'] = addi_muon.time/100 # "time" is actually z in cm
  df.loc[addi_muon.index, 'time'] = np.nan

  # For ID=75/76 we need to offset x/y according to the observation level of the
  # matching 5/6 in the next row
  for pid in [75, 76]:
    keep_idx,keep_vals,drop_idx = _find_paired_data(df, pid)
    _obs_level = keep_vals.obs_level
    birth_muon = df.loc[keep_idx, :]
    df.loc[birth_muon.index, 'obs_level'] = _obs_level.values
    df.loc[birth_muon.index, 'x'] = birth_muon.x - obslev_xy['x'][_obs_level - 1]
    df.loc[birth_muon.index, 'y'] = birth_muon.y - obslev_xy['y'][_obs_level - 1]

  # For ID=95/96: gen number g, muon fate index f: g×10 + f
  fated_muon = df.query("95 <= particle_id <= 96")
  _fate_index = fated_muon.generation % 10
  _generation = (fated_muon.generation - _fate_index) / 10
  df.loc[fated_muon.index, 'fate_index'] = _fate_index.astype(int)
  df.loc[fated_muon.index, 'generation'] = _generation.astype(int)
  # fate index=1: Muon track ends because of decay
  # fate index=2: Muon track ends because of nuclear fatal interaction
  # fate index=3: Muon track ends in update by energy or angular cut

  df['is_muon'] = False
  df.loc[df.query("particle_id in [5, 6, 75, 76, 85, 86, 95, 96]").index, 'is_muon'] = True

  new_col_order = [0, 1, 9, 2, 3, 8, 4, 5, 6, 7, 10, 11, 12, 13]
  df=df[ df.columns[new_col_order] ]

  return df
