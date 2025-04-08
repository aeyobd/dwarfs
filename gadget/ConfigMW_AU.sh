
# Basic code operation

    # LEAN
    SELFGRAVITY
    EXTERNALGRAVITY
    EXTERNALGRAVITY_AGAMA
    AGAMA_UNITS
    
# Gravity options

    TREE_NUM_BEFORE_NODESPLIT=3
    RANDOMIZE_DOMAINCENTER
    
# Softening types and particle types
    NSOFTCLASSES=1
    NTYPES=2

# Floating point accuracy

    POSITIONS_IN_64BIT
    DOUBLEPRECISION=1 # mixed precision

# Group finding

    # FOF

# Miscellaneous code options
    OUTPUT_POTENTIAL
    OUTPUT_ACCELERATION
    EVALPOTENTIAL

# Parallel options, not needed except for niagara
    #NUMBER_OF_MPI_LISTENERS_PER_NODE=2
