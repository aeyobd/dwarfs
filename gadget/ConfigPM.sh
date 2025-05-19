
# Basic code operation

    # LEAN
    SELFGRAVITY
    # EXTERNALGRAVITY
    # EXTERNALGRAVITY_MW
    
# Gravity options

    PMGRID=512
    TREEPM_NOTIMESPLIT # required
    RANDOMIZE_DOMAINCENTER
    
# Softening types and particle types
    NSOFTCLASSES=1
    NTYPES=2

# Floating point accuracy

    POSITIONS_IN_64BIT
    DOUBLEPRECISION=1 

# Group finding

    # FOF

# Miscellaneous code options

    OUTPUT_POTENTIAL
    EVALPOTENTIAL
    OUTPUT_ACCELERATION


# Parallel options, not needed except for niagara
    NUMBER_OF_MPI_LISTENERS_PER_NODE=2
