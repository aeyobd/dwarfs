
# Basic code operation
    SELFGRAVITY
    
# Gravity options
    TREE_NUM_BEFORE_NODESPLIT=3
    RANDOMIZE_DOMAINCENTER

    FMM
    MULTIPOLE_ORDER=4
    
# Softening types and particle types
    NSOFTCLASSES=1
    NTYPES=2

# Floating point accuracy
    POSITIONS_IN_64BIT
    DOUBLEPRECISION=1 

# Miscellaneous code options

    OUTPUT_POTENTIAL
    EVALPOTENTIAL
    OUTPUT_ACCELERATION


# Parallel options, not needed except for niagara
    #NUMBER_OF_MPI_LISTENERS_PER_NODE=2
