@testset "Test loc_to_loc3x for different unit cells of the mapleleaf lattice." begin

    # Define the hexagonal unit cell for the mapleleaf lattice.
    uc_hex = UnitCell(;
        lattice_vecs=[[1.5 * sqrt(3), -0.5], [sqrt(3), 2.0]],
        basis_vecs=[
            [0.0, 0.0],
            [sqrt(3) / 2, -1 / 2],
            [sqrt(3), 0],
            [sqrt(3), 1.0],
            [sqrt(3) / 2, 3 / 2],
            [0.0, 1.0],
        ],
    )

    # Define a triangular unit cell for the mapleleaf lattice.
    uc_tri = UnitCell(;
        lattice_vecs=[[1.5 * sqrt(3), -0.5], [sqrt(3), 2.0]],
        basis_vecs=[
            [sqrt(3), 0],
            [0.0, 0.0],
            [sqrt(3) / 2, 1 / 2],
            [sqrt(3), -1],
            [sqrt(3) / 2, -1 / 2],
            [sqrt(3), 1],
        ],
    )

    # Iterate over the defined unit cells (hexagonal and triangular).
    for uc in (uc_hex, uc_tri)
        # Generate the 3x extended unit cell.
        uc3x = get_unitcell_3x(uc)

        # Test that converting locations from the original unit cell to the 3x extended cell works correctly.
        for i in -10:10, j in -10:10, b in eachindex(uc.basis_vecs)
            # Create a location tuple consisting of a lattice position [i, j] and a basis vector index b.
            loc = ([i, j], b)

            # Convert the location to a real space position using the original unit cell.
            v = loc_to_pos(loc[1], loc[2], uc)

            # Convert the location to the 3x extended unit cell.
            loc3x = CMFT.loc_to_loc3x(loc, uc)

            # Convert the new location to a real space position using the extended unit cell.
            v3x = CMFT.loc_to_pos(loc3x[1], loc3x[2], uc3x)

            # Verify that the positions in the original unit cell and the 3x extended unit cell match.
            @test v ≈ v3x
        end
    end
end
