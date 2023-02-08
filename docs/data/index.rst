Data Products & Specification
================================

.. include:: /global.rst

.. _data-overview:

Overview of the ELUCID Data
----------------------------

The root directory of the ELUCID run has the following layout:

.. code-block:: text 
    :class: all-code-block

    simulation.hdf5
    output/
        |- snapdir_{i}/
            |- snapshot_{i}.{j}.hdf5
        |- groups_{i}/
            |- tab_{i}.{j}.hdf5
        |- subhalos_{i}/
            |- tab_{i}.{j}.hdf5
    postprocessing/
        |- trees/
            |- SubLink/
                |- tree.{j}.hdf5

Where ``i`` is the snapshot number ranging in ``[0, NumSnapshots-1]``. Data in a snapshot 
is spread into ``NumFilesPerSnapshot`` files (i.e., chunks), 
and ``j`` is the file number ranging in ``[0, NumFilesPerSnapshot-1]``. 

Currently there are 101 snapshots (``NumSnapshots=101``) and 2048 files per 
snapshot (``NumFilesPerSnapshot=2048``).

All the data files are organized use HDF5 data format. The details of these files 
are described in the following sections. Below we list a summary of them:

.. table::
    :class: fix-width-table tight-table
    :widths: 30 70

    =================================== ============================================================================================================================
    File or Directory                   Description
    =================================== ============================================================================================================================
    ``simulation.hd5``                  The "all-in-one" interface file. Data in other files are all mapped into this file using HDF5 Virtual Dataset. 
                                        Hence, you may visit all the data without explicitly touching any other file.
    ``output/snapdir_{i}``/             Particle data in the snapshot ``i``. The particles are spread into ``NumFilesPerSnapshot`` separate files according to the 
                                        Peano-Hilbert curves, so that each file stores particles in adjacent cells.
    ``output/groups_{i}``/              Group table in the snapshot ``i``. Each file gives the cell indices and particle indices of a list of groups.
    ``output/subhalos_{i}``/            Subhalo table in the snapshot ``i``. Each file gives the cell indices and particle indices of a list of subhalos.
    ``postprocessing/trees/SubLink/``   The SubLink subhalo merger trees.
    =================================== ============================================================================================================================

.. _data-overview-simulation-hdf5-file:

The simulation.hdf5 File
"""""""""""""""""""""""""""

The ``simulation.hdf5`` file provides a "all-in-one" interface for all the simulation data in this run. 
All the data in this run are mapped into this single file using HDF5 Virtual Dataset. 

It is recommended to only use this file, because the underlying library is responsible for reliably and efficiently 
retreiving the data from the file systems. The users just manipulate the datagroups, datasets and attributs.

The datagroup and dataset layout in ``simulation.hdf5`` is:

.. code-block:: text
    :class: all-code-block

    Header [__attrs__ = BoxSize, HubbleParam, Omega0, OmegaLambda, 
            NumSnapshots, NumFilesPerSnapshot, NumPart_Total, MassTable, 
            GroupTableMaskNumBit, HashTableSize,
            FlagCooling, FlagFeedback, FlagMetals, FlagSFR, FlagStellarAge] /
        |- Redshifts, SnapNums, Times
    Snapshots/{i} [__attrs__ = Redshift, Time, NumGroup_Total] / {j}/
        |- FoFTable [__attrs__ = NumGroup_ThisFile, NumPart_ThisFile] /
            |- GroupLen, GroupOff
            |- GroupParticleIDs
        |- SubfindTable [__attrs__ = NumFilesPerSnapshot, NumGroup_ThisFile, 
                NumGroup_Total, NumPart_ThisFile, NumSub_ThisFile] /
            |- FirstSubOfGroup, NumSubPerGroup
            |- Halo_M_Crit200, Halo_M_Mean200, Halo_M_TopHat200
            |- Halo_R_Crit200, Halo_R_Mean200, Halo_R_TopHat200
            |- SubOff, SubLen, SubParentHalo, SubMostBoundID
            |- SubHalfMass, SubPos, SubSpin, SubVel, SubVelDisp, SubVmax
            |- ParticleIDs
        |- HashTableType1 [__attrs__ = FirstCellID, LastCellID] /
            |- ParticleOffsetsInCell
        |- PartType1 [__attrs__ = NumPart_ThisFile] /
            |- Coordinates, Velocities, ParticleIDs
    Trees/SubLink/{j}/ 
        |- Header [__attrs__ = NumHalos, NumTrees] /
            |- NumHalosInTree, HaloOffsetInTree
        |- Subhalos/
            |- Descendant, FirstProgenitor, NextProgenitor
            |- FirstHaloInFOFgroup, NextHaloInFOFgroup
            |- M_Crit200, M_Mean200, M_TopHat
            |- Pos, Vel
            |- SnapNum, FileNr, SubhaloIndex, Len, MostBoundID
            |- VelDisp, Vmax, SubhalfMass, Spin

Note that: 

- ``i`` is the snapshot number and ``j`` is the file (i.e., chunk) number. 
- ``__attrs__`` means the data object has a set of attached attributes that hold
  the metadata, such as the simulation parameters, number of objects in this 
  data group, etc.

We will demonstrate the data I/O on ``simulation.hdf5`` 
in the following sections using ``h5py``, a Python library for manipulating HDF5 files. 
In addition, we also need ``pandas`` and ``numpy``, as they are extensively used 
in data analysis. To use them, simply import the API by::

    >>> import h5py
    >>> import pandas as pd
    >>> import numpy as np

The C/C++ examples are described in :ref:`data-using-c-cpp`.

.. _data-overview-simulation-parameters:

Simulation Parameters
"""""""""""""""""""""""

The simulation settings and cosmological parameters are:

.. table:: 
    :class: fix-width-table tight--table
    :widths: 35 15 35 15

    =========================================================== ===================================== =========================================================== =================================================
    Parameter                                                   Value                                 Parameter                                                   Value
    =========================================================== ===================================== =========================================================== =================================================
    Number of snapshots :math:`N_{\rm snapshot}`                101                                   Number of files per snapshot :math:`N_{\rm files}`          2048                                  
    Bos size :math:`L_{\rm box}`                                :math:`500.0 h^{-1}\, {\rm Mpc}`      Number of space-filling cells :math:`N_{\rm hash}`          16777216                              
    Number of dark matter particles :math:`N_{\rm dm}`          28991029248                           Dark matter particle mass :math:`M_{\rm dm}`                :math:`0.03087502 \times 10^{10} h^{-1}\,M_\odot`
    Hubble parameter :math:`h`                                  0.72
    :math:`\Omega_{0}`                                          0.258                                 :math:`\Omega_{\Lambda}`                                    0.742
    Group table mask number of bits                             36
    =========================================================== ===================================== =========================================================== =================================================

The redshifts and scale factors for those snapshots are:

.. table::
    :class: fix-width-table tight-table
    :widths: 10 10 13 10 10 13 10 10 13
    
    ============== ========= ============= ========== ========= ========= ============ =========  ==========
    Snapshot       :math:`z` :math:`a`     Snapshot   :math:`z` :math:`a` Snapshot     :math:`z`  :math:`a`                                                    
    ============== ========= ============= ========== ========= ========= ============ =========  ==========
    0              18.409561 0.051521      34         6.009231  0.142669  68           1.531159   0.395076
    1              17.837003 0.053087      35         5.802351  0.147008  69           1.456453   0.407091
    2              17.280867 0.054702      36         5.601575  0.151479  70           1.383961   0.419470
    3              16.741506 0.056365      37         5.406766  0.156085  71           1.313599   0.432227
    4              16.217927 0.058079      38         5.217668  0.160832  72           1.245319   0.445371
    5              15.709555 0.059846      39         5.034165  0.165723  73           1.179053   0.458915
    6              15.216655 0.061665      40         4.856104  0.170762  74           1.114742   0.472871
    7              14.737870 0.063541      41         4.683271  0.175955  75           1.052330   0.487251
    8              14.273472 0.065473      42         4.515537  0.181306  76           0.991758   0.502069
    9              13.822720 0.067464      43         4.352746  0.186820  77           0.932976   0.517337
    10             13.385178 0.069516      44         4.194778  0.192501  78           0.875930   0.533069
    11             12.960631 0.071630      45         4.041466  0.198355  79           0.820565   0.549280
    12             12.548667 0.073808      46         3.892679  0.204387  80           0.766834   0.565984
    13             12.148725 0.076053      47         3.748270  0.210603  81           0.714689   0.583196
    14             11.760799 0.078365      48         3.608146  0.217007  82           0.664085   0.600931
    15             11.384054 0.080749      49         3.472132  0.223607  83           0.614971   0.619206
    16             11.018653 0.083204      50         3.340146  0.230407  84           0.567310   0.638036
    17             10.663984 0.085734      51         3.212051  0.237414  85           0.521054   0.657439
    18             10.319644 0.088342      52         3.087756  0.244633  86           0.476161   0.677433
    19             9.985631  0.091028      53         2.967105  0.252073  87           0.432595   0.698034
    20             9.661436  0.093796      54         2.850019  0.259739  88           0.390316   0.719261
    21             9.346719  0.096649      55         2.736404  0.267637  89           0.349284   0.741134
    22             9.041370  0.099588      56         2.626131  0.275776  90           0.309461   0.763673
    23             8.745069  0.102616      57         2.519107  0.284163  91           0.270816   0.786896
    24             8.457427  0.105737      58         2.415254  0.292804  92           0.233310   0.810826
    25             8.178270  0.108953      59         2.314452  0.301709  93           0.196911   0.835484
    26             7.907416  0.112266      60         2.216634  0.310884  94           0.161586   0.860892
    27             7.644537  0.115680      61         2.121703  0.320338  95           0.127304   0.887072
    28             7.389403  0.119198      62         2.029569  0.330080  96           0.094034   0.914048
    29             7.141798  0.122823      63         1.940156  0.340118  97           0.061746   0.941845
    30             6.901516  0.126558      64         1.853385  0.350461  98           0.030411   0.970487
    31             6.668300  0.130407      65         1.769170  0.361119  99           0.000000   1.000000
    32             6.442027  0.134372      66         1.687450  0.372100  100          0.000000   1.000000
    33             6.222355  0.138459      67         1.608133  0.383416      
    ============== ========= ============= ========== ========= ========= ============ =========  ==========


All these global parameters are stored in a header group of the simulation file. To get them, open the simulation file 
from Python by creating a ``h5py.File`` object::

    >>> f = h5py.File('simulation.hdf5', 'r')

where the two positional parameters are the file path and the opening/access mode (``r`` for read-only on a existing file).

The names of data groups under the file are::

    >>> f.keys()
    <KeysViewHDF5 ['Header', 'Snapshots', 'Trees']>

To get the simulation parameters, visit the ``Header`` group::

    >>> header = f['Header']
    >>> header.keys()
    <KeysViewHDF5 ['Redshifts', 'SnapNums', 'Times']>

Three datasets are contained in the header. For example, to get the redshifts of all the snapshots, run 

    >>> header['Redshifts'][()]
    array([ 1.8409561e+01,  1.7837003e+01,  1.7280867e+01,  1.6741506e+01, ...])


Other scalar or small-vector parameters are store in the attributes of the header to enable more 
compact storage and more efficient loading::

    >>> header.attrs.keys()
    <KeysViewHDF5 ['BoxSize', 'FlagCooling', 'FlagFeedback', 'FlagMetals', 'FlagSFR', 'FlagStellarAge', 'GroupTableMaskNumBit', 'HashTableSize', 'HubbleParam', 'MassTable', 'NumFilesPerSnapshot', 'NumPart_Total', 'NumSnapshots', 'Omega0', 'OmegaLambda']>

To get those parameters, for example, the Hubble parameter :math:`h` and dark energy density parameter :math:`\Omega_{\Lambda, 0}`, run::

    >>> header.attrs['HubbleParam'], header.attrs['OmegaLambda']
    (0.72, 0.742)

To get the total number of dark matter particles at every snapshot and the particle mass, run

    >>> header.attrs['NumPart_Total'][1], header.attrs['MassTable'][1]
    (28991029248, 0.030875024131732435)

Here the index ``1`` refers to the dark matter particles. Other particles are not available in a dark-matter-only simulation.

.. _data-particle-catalog:

Particle Catalog
-----------------

The three available fields for dark matter particles are 

.. table::
    :class: fix-width-table tight-table
    :widths: 15 15 15 55

    ======================== ============================== ============================= =================================================================================================================================
    Field                    Data type and dimensions       Unit                          Description
    ======================== ============================== ============================= =================================================================================================================================
    ``Coordinates``           ``float32 (N,3)``             :math:`h^{-1}{\rm Mpc}`       Comoving spatial coordinate in periodic box. |br| The value are not bound to ``[0, BoxSize)`` due to periodic boundary condition.
    ``Velocities``            ``float32 (N,3)``             :math:`\sqrt{a}\, {\rm km/s}` Spatial velocity. |br| Multiply this value by :math:`\sqrt{a}` to obtain the peculiar velocity.
    ``ParticleIDs``           ``int64 (N,)``                \-                            Unique ID of this particle. |br| Constant across snapshots. 
    ======================== ============================== ============================= =================================================================================================================================


Particles in the snapshot ``i`` are stored in the data group ``Snapshots/{i}`` as ``NumFilesPerSnapshot`` chunks separately. For 
example, to get the the first chunk of particles in the snapshot ``76``, run::

    >>> parts = f['Snapshots/76/0/PartType1']
    >>> parts.keys()
    <KeysViewHDF5 ['Coordinates', 'ParticleIDs', 'Velocities']>

The number of dark matter particles in this file is::
    
    >>> parts.attrs['NumPart_ThisFile']
    13158492

To load the position of, for example, the first 10 particles in this chunk, run::

    >>> parts['Coordinates'][:10]
    array([[0.2370728 , 0.4263347 , 0.282894  ],
       [0.2030037 , 0.29562455, 0.30211562],
       [0.18131527, 0.06854798, 0.25142768],
       [0.16527686, 0.11370011, 0.7375622 ],
       [0.14694212, 0.2808326 , 0.7435081 ],
       [0.13858652, 0.41486353, 0.7316461 ],
       [0.66691387, 0.23588769, 0.6613449 ],
       [0.70934445, 0.06762647, 0.7157361 ],
       [0.68830323, 0.00258765, 0.712928  ],
       [0.6061265 , 0.39591527, 0.58704   ]], dtype=float32)

.. _data-group-table:

Group Table
--------------

The group table stores the output of FoF algorithm, i.e., how many particles are linked to each FoF group and the IDs of those particles.
The combination of group table and particle catalog allows extracting particles associated with each FoF group.

The available fields for groups are:

.. table::
    :class: fix-width-table tight-table
    :widths: 20 20 60

    =================== =============================== ===============================================================================================
    Field               Data type and dimensions        Description
    =================== =============================== ===============================================================================================
    GroupLen            ``int32, (Ngroup, )``           Number of particles associated with each FoF group.
    GroupOff            ``int32, (Ngroup, )``           The offset (i.e., index to ``GroupParticleIDs``) of the first particle in each FoF group. 
    GroupParticleIDs    ``int64, (Nparticle, )``        The cell ID and particle ID of each particle.
    =================== =============================== ===============================================================================================

Group table in the snapshot ``i`` is stored in the data group ``Snapshots/{i}`` as ``NumFilesPerSnapshot`` chunks separately. 
Each chunk stores a distinct set of groups, i.e., the information of one group cannot be across two chunks.
However, the chunks in the group table do not one-to-one match the chunks in the particle data, 
i.e., the particles associated with groups 
in the chunk ``j`` of group table may be distributed in multiple chunks in the particle data.

For example, to get the first chunk of the group table in the snapshot ``76``, run::

    >>> groups = f['Snapshots/76/0/FoFTable']
    >>> groups.keys()
    <KeysViewHDF5 ['GroupLen', 'GroupOff', 'GroupParticleIDs']>

The number of groups in this chunk and the total number of particles **associated with them** are::

    >>> groups.attrs['NumGroup_ThisFile'], groups.attrs['NumPart_ThisFile']
    (24003, 4708732)

The number of particles in each of the first 10 groups are::

    >>> groups['GroupLen'][:10]
    array([305341, 132169,  95756,  87657,  79336,  60440,  58894,  57956,
        56591,  48662], dtype=int32)

To get the particle indices, for example, in the 10th ``(index=9)`` group, combind the ``GroupOff`` and ``GroupLen`` fields as indexing 
range into ``GroupParticleIDs``::

    >>> offset, length = groups['GroupOff'][9], groups['GroupLen'][9]
    >>> pids = groups['GroupParticleIDs'][offset:offset+length]
    >>> pids
    array([251994538396849, 252063059717285, 252063059717286, ...,
        252819417542811, 252819426979993, 252819426979994])

Note that the ``pids`` are bit-wise ``OR`` ed from cell IDs and real particle IDs 
(to save disk and memory storage). 
Specifically, the cell ID of each particle is logically left shifted ``GroupTableMaskNumBit``
bits, and then ``OR`` ed with the particle ID. 
The number of bits shifted is stored in the simulation header::

    >>> mask_nbit = f['Header'].attrs['GroupTableMaskNumBit']
    >>> mask_nbit
    36

To retrieve the cell IDs, run::

    >>> cids = pids >> mask_nbit
    >>> cids
    array([3667, 3668, 3668, ..., 3679, 3679, 3679])

To retrieve the (real) particle IDs, run::

    >>> mask = (1 << mask_nbit) - 1
    >>> pids = pids & mask
    >>> pids
    array([217205937,  19049637,  19049638, ..., 462631067, 472068249, 472068250])


The cell IDs and particle IDs, if combined with the ``HashTableType1`` and 
``PartType1`` data groups, allow us to load the particles belonging to this 
group. Note that the particles of a group can be separately stored 
in different chunks (files). To work out, we can proceed as follows:

- Build a mapping from cell index to chunk (file) index. From that, we 
  map ``cids`` to ``fids``, the indices of files that contain these cells.
- We open each file in ``fids``. For each file, we load the data in each cell
  in ``cids`` that are contained in this file.

In Python, we may implement the loading process as follows. 

The ``cid_to_fid`` mapping can be constructed by visiting the hash table in each 
file, using the staring and ending cell indices in that file::

    # Build an array 'cid_to_fid' that maps cell index to file index
    header = f['Header']
    n_files, n_cells = header.attrs['NumFilesPerSnapshot'], header.attrs['HashTableSize']
    cid_to_fid = np.empty(n_cells, dtype=int)

    for j in range(n_files):
        hash_tab = f[f'Snapshots/76/{j}/HashTableType1']
        first_cid, last_cid = hash_tab.attrs['FirstCellID'], hash_tab.attrs['LastCellID']
        cid_to_fid[first_cid:last_cid+1] = j


Then, we find the file indices for all target cells::

    cids = sorted(set(cids))
    fids = cid_to_fid[cids]
    pids = pd.Index(pids)
    cids, fids, pids

Here, the ``set`` class is used to drop duplicated cells. The sorting operation 
ensures the monotonical ordering of the resulted ``fids``, from which we can easily 
visit the files one by one and avoid repeated close-and-reopen.
The pandas index type, ``pd.Index``, will be used later for fast intersection 
operation. The output is:

.. code-block:: text
    :class: all-code-block

    [3667, 3668, 3671, 3672, 3675, 3676, 3679],
    array([0, 0, 0, 0, 0, 0, 0]),
    Int64Index([217205937,  19049637,  19049638,  19052707, ...],
        dtype='int64', length=48662)  

We see that the particles of the group are contained in 7 cells, and that 
all these cells are contained in the first file (indexed 0). We then iterate over
each required file, and in each file, we iterate over all required cells::

    cur_fid = -1
    coords = []
    for fid, cid in zip(fids, cids):
        if cur_fid != fid:
            cur_fid = fid
            
            # We only open the dataset handles, but do not load the data immediately
            # for the sake of performance.
            parts = f[f'Snapshots/76/{fid}/PartType1']
            coords_in_file = parts['Coordinates']
            pids_in_file = parts['ParticleIDs']
            n_ps = coords_in_file.shape[0]
            
            hash_tab = f[f'Snapshots/76/{fid}/HashTableType1']
            first_cid, last_cid = hash_tab.attrs['FirstCellID'], hash_tab.attrs['LastCellID']
            poffs = np.hstack((hash_tab['ParticleOffsetsInCell'][()], n_ps)) 
            
            print(f'File {fid} loaded (n_particles={n_ps}, cell_id range=[{first_cid}, {last_cid}])')
            
        # For the cell indexed 'cid', we find the offset range of its particles, 
        # in this file indexed 'fid'.
        # We then load the information of particles in that cell.
        first_poff, last_poff = poffs[cid-first_cid], poffs[cid-first_cid+1]
        coords_in_cell = coords_in_file[first_poff:last_poff]
        pids_in_cell = pids_in_file[first_poff:last_poff]
        
        # Only a subset of particles in the cell belongs to the target group. 
        # We use the Index instance to find them and append them to the output.
        mask = pids.get_indexer(pids_in_cell) >= 0
        coords_in_target = coords_in_cell[mask]
        coords.append(coords_in_target)
        print(f'In cell {cid}, {len(coords_in_target)} particles found')
        
    coords = np.concatenate(coords)
    assert(coords.shape == (pids.size, 3))      # Make sure that all particles were found.

The ``coords`` array will eventually contains the positions of particles in 
the target group. The output is:

.. code-block:: text
    :class: all-code-block

    File 0 loaded (n_particles=13158492, cell_id range=[0, 8191])
    In cell 3667, 1 particles found
    In cell 3668, 8779 particles found
    In cell 3671, 2292 particles found
    In cell 3672, 9356 particles found
    In cell 3675, 634 particles found
    In cell 3676, 682 particles found
    In cell 3679, 26918 particles found

Altough we have 8192 cells in this file, we loaded only 7 cells out of them. 
Given a small sample of target groups/subhalos/particles, this strategy largely 
reduces the amount of data to be loaded. Depending on the target sample 
in your application, the details of the loading algorithm should be tuned 
if the IO performance is a issue.

.. _data-subhalo-table:

Subhalo Table
--------------

The subhalo table stores the output of the Subfind algorithm, i.e., the subhalo memberships
of each group, and the particle memberships of each subhalo.
In addition, several group properties, such as virial radius and mass, and subhalo 
properties, such as position and velocity, are also provided.

The combination of subhalo table and particle catalog allows extracting particles 
associated with each Subfind subhalo.

Available datasets for the subhalo table are:

.. table::
    :class: fix-width-table tight-table
    :widths: 20 20 60

    =================== =============================== ===============================================================================================
    Field               Data type and dimensions        Description
    =================== =============================== ===============================================================================================
    FirstSubOfGroup     ``int32, (Ngroup, )``           The index (into subhalo datasets) of the first subhalo of each group.
    NumSubPerGroup      ``int32, (Ngroup, )``           Number of subhalos of each group. 
    SubOff              ``int32, (Nsubhalo, )``         The index into ``ParticleIDs`` of the first particle of each subhalo. 
    SubLen              ``int32, (Nsubhalo, )``         Number of particles associated with each subhalo.
    SubParentHalo       ``int32, (Nsubhalo, )``         The index (into group datasets) of the host group of each subhalo.
    SubMostBoundID      ``int64, (Nsubhalo, )``         The particle ID of the most bound particle of each subhalo (no bit-mask is needed).
    ParticleIDs         ``int64, (Nparticle, )``        The cell ID and particle ID of each particle (combined with bit-mask).
    =================== =============================== ===============================================================================================

Other datasets not listed here are ``Halo_M_Crit200``, ``Halo_M_Mean200``, ``Halo_M_TopHat200``, ``Halo_R_Crit200``, 
``Halo_R_Mean200``, ``Halo_R_TopHat200``, ``SubHalfMass``, ``SubPos``, ``SubSpin``, 
``SubVel``, ``SubVelDisp``, and ``SubVmax``. See the merger tree catalog below for their detailed meanings.

Similar to dealing with group tables, for example, to load the first chunk of 
the subhalo table at snapshot ``99``, we write::

    >>> subs = f['Snapshots/99/0/SubfindTable']

The meta info, such as numbers of groups, subhalos and particles in this chunk, 
can be retrieved from the attributes of the data group::

    >>> attrs = subs.attrs
    >>> attrs['NumGroup_ThisFile'], attrs['NumSub_ThisFile'], attrs['NumPart_ThisFile']
    (24873, 25500, 7440527)

To find information of all subhalos in the 9-th group, for example, their maximal 
circular velocities, run::

    >>> grp_id = 9
    >>> first_sub_id = subs['FirstSubOfGroup'][grp_id]
    >>> n_subs = subs['NumSubPerGroup'][grp_id]
    >>> subs['SubVmax'][first_sub_id:first_sub_id+n_subs]    
    array([291.05557 , 353.88458 , 231.2501  , 186.4264  , 164.374   ,
        133.07835 , 117.472244, 141.52196 , 113.723656, 139.05486 , ...])

To locate all particles in the first subhalo of this group, write::

    >>> offset, length = subs['SubOff'][first_sub_id], subs['SubLen'][first_sub_id]
    >>> subs['ParticleIDs'][offset:offset+length]
    array([75180079984635, 75180098859016, 75180061119494, ...,
       75179881763842, 75179957374951, 75179872320523])

Note that the printed "indices" are combination of cell IDs and particles 
IDs. To extract the actual particle IDs, use the same masking method as 
described above for the :ref:`data-group-table`.

.. _data-merger-trees:

Merger Trees
--------------

Subhalos in the ELUCID run are linked by the SubLink algorithm across different redshifts. 
The linked subhalos are represented as a tree called subhalo merger tree which gives the 
assembly histories of the subhalos on the tree.

We use the three-pointer tree representation to store the merger tree, which is described
in any algorithm textbook. Specifically:

- The ``FirstProgenitor`` of a subhalo 
  ``h`` is the index of the most massive progenitor subhalo of ``h``. 
- The ``NextProgenitor``
  of ``h`` is the index of the next most massive progenitor which shares the same descendant 
  as ``h``. 
- The ``Descendant`` of ``h`` is the index of the descendant subhalo of ``h``.

For convenience, two additional indices are provided as linkages of subhalos in the 
same FoF group:

- The ``FirstHaloInFOFgroup`` of ``h`` is the index of the 
  central subhalo (i.e., ``?????`` most massive, as 
  calculated by sum of linked dark matter particles) in the same FoF group as ``h``.
  Hence, a central subhalo has ``FirstHaloInFOFgroup`` pointing to itself. 
- The ``NextHaloInFOFgroup`` of ``h`` is the index of the next most massive subhalo 
  in the same FoF group as ``h``.

Subhalos connected by any of the above five indices are considered in the same tree,
and they are stored contiguously.
All the five indices are ``0-based`` and started
from the first subhalo in each tree. All of them, except ``FirstHaloInFOFgroup``,
can be ``-1``, indicating none of such link exists. Note that the ``FirstProgenitor``
of ``Descendant`` chain may skip at most one snapshot because of those unresolved 
subhalos.

The following list gives detailed description of the fields for subhalos in the tree:

.. table::
    :class: fix-width-table tight-table
    :widths: 15 15 15 55

    ======================== ============================== ================================================= =================================================================================================================================
    Field                    Data type and dimensions       Unit                                              Description
    ======================== ============================== ================================================= =================================================================================================================================
    ``FirstProgenitor``      ``int32, (N,)``                \-                                                Index of the most massive progenitor subhalo.
    ``NextProgenitor``       ``int32, (N,)``                \-                                                Index of the next most massive progenitor subhalo that shares the same descendant.
    ``Descendant``           ``int32, (N,)``                \-                                                Index of the descendant subhalos.
    ``FirstHaloInFOFgroup``  ``int32, (N,)``                \-                                                Index of the most massive subhalo in the same FoF group.
    ``NextHaloInFOFgroup``   ``int32, (N,)``                \-                                                Index of the next most massive subhalo in the same FoF group.
    ``Pos``                  ``float32, (N, 3)``            :math:`h^{-1}\,{\rm Mpc}`                         Comoving spatial position in the periodic box, defined as the position of its most bound particle.
    ``Vel``                  ``float32, (N, 3)``            :math:`\rm km/s`                                  Peculiar velocity of the subhalo, defined as the averaged velocities of all particles linked to it.
    ``M_Mean200``            ``float32, (N,)``              :math:`10^{10}\,h^{-1}M_\odot`                    Total mass of the FoF group enclosed in a sphere whose mean density is 200 times the mean density of the Universe of that time.
                                                                                                              This field is only significant for a central subhalo and set to ``0`` for a satellite subhalo.
    ``M_Crit200``            ``float32, (N,)``              :math:`10^{10}\,h^{-1}M_\odot`                    Total mass of the FoF group enclosed in a sphere whose mean density is 200 times the critical density of the Universe of that time.
                                                                                                              This field is only significant for a central subhalo and set to ``0`` for a satellite subhalo.
    ``M_TopHat``             ``float32, (N,)``              :math:`10^{10}\,h^{-1}M_\odot`                    Total mass of the FoF group enclosed in a sphere whose mean density is :math:`\Delta_c` times the critical density of the 
                                                                                                              Universe of that time, according to the spherical collapse model of Bryan+1998.
                                                                                                              This field is only significant for a central subhalo and set to ``0`` for a satellite subhalo.
    ``Vmax``                 ``float32, (N,)``              :math:`\rm km/s`                                  Maximal peculiar velocity of the spherical-averaged rotation curve.
    ``Len``                  ``int32, (N,)``                \-                                                Number of particles linked to this subhalo.
    ``Spin``                 ``float32, (N, 3)``            :math:`h^{-1}{\rm Mpc\,km/s}`                     Total spin per axis, computed as the average of relative position times relative velocity of all linked particles.
    ``VelDisp``              ``float32, (N,)``              :math:`\rm km/s`                                  1-D velocity dispersion of all linked particles (3-D dispersion divided by :math:`\sqrt{3}`).
    ``SubhalfMass``          ``float32, (N,)``              \-                                                ``?????????????``
    ``SnapNum``              ``int32, (N,)``                \-                                                Snapshot number of this subhalo.
    ``FileNr``               ``int32, (N,)``                \-                                                ``????``
    ``SubhaloIndex``         ``int32, (N,)``                \-                                                ``????``
    ``MostBoundID``          ``int64, (N,)``                \-                                                The particle ID of the most bound particle linked with the subhalos.
    ======================== ============================== ================================================= =================================================================================================================================


SubLink trees are stored in the data group ``Trees/SubLink/{i}`` as ``NumFilesPerSnapshot`` chunks separately. 
Each chunk stores a distinct set of trees, i.e., a forest. One tree cannot be across two chunks. 
For each tree, its subhalos are stored contiguously.

For example, the first chunk can be visited by::

    >>> forest = f['Trees/SubLink/0']
    >>> header, subhs = forest['Header'], forest['Subhalos']
    >>> header, subhs
    (<HDF5 group "/Trees/SubLink/0/Header" (2 members)>,
    <HDF5 group "/Trees/SubLink/0/Subhalos" (19 members)>)

Here, the metainfo are subhalos are contained in two data groups, respectively.

The number of trees and number of subhalos in this chunk is::

    >>> header.attrs['NumTrees'], header.attrs['NumHalos']
    (18326, 1715555)

To get the ``FirstProgenitor`` of subhalos in the 10th tree (i.e., ``index=9`` ),
first load the offset of the first subhalo in the tree, and number of subhalos 
in the tree, by::

    >>> off, nhalos = header['HaloOffsetInTree'][9], header['NumHalosInTree'][9]
    >>> off, nhalos
    (590539, 14210)

Then, load the required subhalo field::

    >>> fpros = subhs['FirstProgenitor'][off:off+nhalos]
    >>> fpros
    array([    1,     2,     3, ..., 14208, 14209,    -1], dtype=int32)

The tophat mass of the ``FirstProgenitor`` of the first subhalo is::

    >>> subhs['M_TopHat'][off + fpros[0]]
    1169.5768

.. _data-in-shell:

Navigate the Data in Shell
-----------------------------------

HDF5 library officially provides a set of command line tools. For example, ``h5ls`` 
command enables a quick navigation of data in any HDF5 file.

For example, to list data groups or datasets in a file, run:

.. code-block:: bash

    $ h5ls simulation.hdf5
    Header                   Group
    Snapshots                Group
    Trees                    Group

To look at all data groups and datasets below a certain level, use the recursive option ``-r``, 
and a ``/`` separarated path followed by the file name 
to specify the target data group:

.. code-block:: bash

    $ h5ls -r simulation.hdf5/Trees/SubLink/0
    /Header                  Group
    /Header/HaloOffsetInTree Dataset {18326}
    /Header/NumHalosInTree   Dataset {18326}
    /Subhalos                Group
    /Subhalos/Descendant     Dataset {1715555}
    /Subhalos/FileNr         Dataset {1715555}
    /Subhalos/FirstHaloInFOFgroup Dataset {1715555}
    /Subhalos/FirstProgenitor Dataset {1715555}
    /Subhalos/Len            Dataset {1715555}
    /Subhalos/M_Crit200      Dataset {1715555}
    /Subhalos/M_Mean200      Dataset {1715555}
    /Subhalos/M_TopHat       Dataset {1715555}
    /Subhalos/MostBoundID    Dataset {1715555}
    /Subhalos/NextHaloInFOFgroup Dataset {1715555}
    /Subhalos/NextProgenitor Dataset {1715555}
    /Subhalos/Pos            Dataset {1715555, 3}
    /Subhalos/SnapNum        Dataset {1715555}
    /Subhalos/Spin           Dataset {1715555, 3}
    /Subhalos/SubhalfMass    Dataset {1715555}
    /Subhalos/SubhaloIndex   Dataset {1715555}
    /Subhalos/Vel            Dataset {1715555, 3}
    /Subhalos/VelDisp        Dataset {1715555}
    /Subhalos/Vmax           Dataset {1715555}

Here, the metainfo of all the datasets in the first chunk of merger tree are printed.

To see the detail contents of datasets or attributes, use ``h5dump`` command. For example,
the ``-a`` option allows printing the attribute detail:

.. code-block:: bash

    $ h5dump -a /Header/HubbleParam simulation.hdf5
    HDF5 "simulation.hdf5" {
        ATTRIBUTE "HubbleParam" {
        DATATYPE  H5T_IEEE_F64LE
        DATASPACE  SCALAR
        DATA {
            (0): 0.72
        }
    }
    }

Here, we know that the Hubble parameter is ``0.72``, stored as a double-precision floating-point scalar value.

To see a dataset, use ``-d`` option followed by the dataset path, ``-s`` to specify the starting indices 
of the dataset array at every dimension, ``-c`` to specify the number of elements to print at every dimension:

.. code-block:: bash

    $ h5dump -d /Snapshots/76/0/PartType1/Coordinates -s '0,0' -c '5,3' simulation.hdf5
    HDF5 "simulation.hdf5" {
    DATASET "/Snapshots/76/0/PartType1/Coordinates" {
       DATATYPE  H5T_IEEE_F32LE
       DATASPACE  SIMPLE { ( 13158492, 3 ) / ( 13158492, 3 ) }
       SUBSET {
          START ( 0, 0 );
          STRIDE ( 1, 1 );
          COUNT ( 5, 3 );
          BLOCK ( 1, 1 );
          DATA {
          (0,0): 0.237073, 0.426335, 0.282894,
          (1,0): 0.203004, 0.295625, 0.302116,
          (2,0): 0.181315, 0.068548, 0.251428,
          (3,0): 0.165277, 0.1137, 0.737562,
          (4,0): 0.146942, 0.280833, 0.743508
          }
       }
    }
    }

Here, we print the spatial positions of the first five dark matter particles in the first chunk of snapshot ``76``.

.. _data-using-c-cpp:

Using ELUCID Data with C/C++ 
----------------------------------

The HDF5 official library provides a set of C-API to manipulate HDF5 files. However, the HPC toolkit set 
`HIPP <https://github.com/ChenYangyao/hipp>`_ provides more elegant ``C++`` API without significant overhead. 
In this section, we demonstrate loading ELUCID dataset with the IO module of HIPP. 

The source files of the examples can be downloaded:

.. table::
    :class: fix-width-table tight-table
    :widths: 30 70

    =================================================================================== ===========================
    Source file                                                                         Description
    =================================================================================== ===========================
    :download:`load-particles.cpp <../examples/cpp/load-particles.cpp>`                 Load particles.
    :download:`load-merger-tree.cpp <../examples/cpp/load-merger-tree.cpp>`             Load merger trees.
    :download:`Makefile <../examples/cpp/Makefile>`                                     The makefile for compiling
    =================================================================================== ===========================

We assume the following header inclusion and namespace declaration for clarity:

.. code-block:: cpp

    #include <hippio.h>

    using namespace HIPP;
    using namespace HIPP::IO;

Loading Particle Data
""""""""""""""""""""""""""

To open a HDF5 file, declare a ``H5File`` object with the file name and opening/access mode passed. 
The ``open_group()`` method opens a data group under a file or a parent data group. For example, 
the header can be opened by:

.. code-block:: cpp

    H5File f("simulation.hdf5", "r");
    H5Group header = f.open_group("Header");

To get access attributes, use ``open_attr()`` method on HDF5 objects (file, group, or dataset), and 
call ``read()`` on the returned attribute object. For example, the cosmological parameters can be obtained 
from the attributes of the header:

.. code-block:: cpp

    double h, Omg0;
    header.open_attr("HubbleParam").read(&h);
    header.open_attr("Omega0").read(&Omg0);
    pout << "hubble parameter = ", h, ", matter density = ", Omg0, endl;

Output:

.. code-block:: text
    
    hubble parameter = 0.72, matter density = 0.258


To load a dataset, call ``open_dataset()`` on a data group, and call ``read()`` on
the returned dataset object. For example, to load redshifts in all snapshots, use:

.. code-block:: cpp

    vector<double> zs;
    header.open_dataset("Redshifts").read(zs);
    pout << "Redshifts = {", zs, "}", endl;

Note that the ``read()`` call allows passing a vector whose size is automatically adapted
or passing a pointer to a scalar object with pre-allocated memory.

Output:

.. code-block:: text

    Redshifts = {18.4096,17.837,17.2809,16.7415, ...}

To load particle data, for example, open the first chunk of snapshot ``76``:

.. code-block:: cpp

    auto parts = f.open_group("Snapshots/76/0/PartType1");
    vector<string> keys = parts.keys();
    pout << "Names of fields: {", keys, "}\n";

Here, we have queried the names of the datasets by ``key()`` on the data group. The output is:

.. code-block:: text

    Names of fields: {Coordinates,ParticleIDs,Velocities}
    
To get the number of particles, visit its attribute:

.. code-block:: cpp

    long long n_parts;
    parts.open_attr("NumPart_ThisFile").read(&n_parts);
    pout << "No. of particles = ", n_parts, endl;

Output:

.. code-block:: text

    No. of particles = 13158492

To load coordinates of all particles in this chunk, open the dataset ``Coordinates``
and read its content into a preallocated memory buffer:

.. code-block:: cpp

    vector<std::array<float, 3>> xs(n_parts);
    parts.open_dataset("Coordinates").read(&xs[0][0]);
    pout << "Positions of the first five particles: \n";
    for(int i=0; i<5; ++i)
        pout << "  {", xs[i], "}\n";

Output:

.. code-block:: text

    Positions of the first five particles:
    {0.237073,0.426335,0.282894}
    {0.203004,0.295625,0.302116}
    {0.181315,0.068548,0.251428}
    {0.165277,0.1137,0.737562}
    {0.146942,0.280833,0.743508}


Loading Merger Trees 
""""""""""""""""""""""

Merger trees can be loaded using the same way as that of the particles. However,
we usually need more than one field. To tackle this, we declare a structure 
type to hold a subhalo:

.. code-block:: cpp

    struct Subhalo {
        float mass;
        float x[3];
        int fpro_id, cent_id;
        int snap;
    };


We open the first chunk of the tree data by:

.. code-block:: cpp

    H5File f("simulation.hdf5", "r");
    auto g = f.open_group("Trees/SubLink/0");

To load multiple fields into an array of a structure type, declare a ``H5XTable``
object, with the structure type as its template parameter, and 
the list of dataset names and member pointers as arguments to the constructor. 
Then, call ``read()`` on the object by passing a target data group. A vector 
of structure objects is returned:

.. code-block:: cpp

    vector<Subhalo> forest = H5XTable<Subhalo>(
        "M_TopHat", &Subhalo::mass,                 "Pos", &Subhalo::x,
        "FirstProgenitor", &Subhalo::fpro_id,       "FirstHaloInFOFgroup", &Subhalo::cent_id,
        "SnapNum", &Subhalo::snap
    ).read(g.open_group("Subhalos"));

Note that we have loaded all trees (i.e., a forest) from the chunk.

To load the tree offsets and the numbers of subhalos in trees, visit the datasets 
in the header:

.. code-block:: cpp

    vector<int> offs, lens;
    g.open_dataset("Header/HaloOffsetInTree").read(offs);
    g.open_dataset("Header/NumHalosInTree").read(lens);

    int nhalos = forest.size(), ntrees = offs.size();
    pout << "No. of subhalos/trees in this chunk = ", nhalos, '/', ntrees, 
        endl;

Output:

.. code-block:: text

    No. of subhalos/trees in this chunk = 1715555/18326

The index range of a tree is determined by the offset (of its first subhalo) and
the number of subhalos. For example, to visit the first tree (``index=0``), simply
get the pointer to its first subhalo by

.. code-block:: cpp

    int off = offs[0], len = lens[0];
    const auto *tree = &forest[off];

    
As an application, we locate a central, ``z=0``, massive subhalo in the tree
by iterating over all subhalos in it:

.. code-block:: cpp
    
    int root_id = -1;
    for(int i=0; i<len; ++i){
        auto &h = tree[i];
        if( h.snap == 100 && h.cent_id == i && h.mass > 1.0e2 ) {
            root_id = i; break;
        }
    }
    pout << "A root index at z=0 : ", root_id, endl;

Output:

.. code-block:: text

    A root index at z=0 : 0
    
Then, we (roughly) estimate the half-mass formation time by following the first-progenitor
chain and stopping until the mass attribute going below half of the root mass:

.. code-block:: cpp

    int pro_id = root_id;
    while( pro_id >= 0 ) {
        auto &pro = tree[pro_id];
        if( pro.mass < tree[root_id].mass * 0.5 ) {
            pout << "Half mass formation snapshot: ", pro.snap, endl;
            break;
        }
        pro_id = pro.fpro_id;
    }

Output:

.. code-block:: text
    
    Half mass formation snapshot: 89


