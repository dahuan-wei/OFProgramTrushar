    scalar CoNum = -GREAT;

    forAll(fluidRegions, regionI)
    {
        CoNum = max
        (
            icoCourantNo
            (
                fluidRegions[regionI],
                runTime,
                phiFluid[regionI]
            ),
            CoNum
        );
    }
/*
    forAll (porousFluidRegions, porousI)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                porousFluidRegions[porousI],
                runTime,
                rhoPorous[porousI],
                phiPorous[porousI]
            ),
            CoNum
        );
    }
*/
