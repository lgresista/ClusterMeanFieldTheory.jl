#digits function for base 2
digits2(i::Integer) = digits(i; base=2, pad=8)

for _ in 0:9
    n = rand(UInt8)
    v = digits2(n)

    # setbit
    for i in 0:7
        np = CMFT.setbit(n, i)

        vp = copy(v)
        vp[i + 1] = 1

        @test digits2(np) == vp
    end

    # unsetbit
    for i in 0:7
        np = CMFT.unsetbit(n, i)

        vp = copy(v)
        vp[i + 1] = 0

        @test digits2(np) == vp
    end

    # flipbit
    for i in 0:7
        np = CMFT.flipbit(n, i)

        vp = copy(v)

        if vp[i + 1] == 0
            vp[i + 1] = 1
        else
            vp[i + 1] = 0
        end

        @test digits2(np) == vp
    end

    # flipbits
    for i in 0:7
        for j in 0:7
            np = CMFT.flipbits(n, i, j)

            vp = copy(v)
            if vp[i + 1] == 0
                vp[i + 1] = 1
            else
                vp[i + 1] = 0
            end

            if vp[j + 1] == 0
                vp[j + 1] = 1
            else
                vp[j + 1] = 0
            end

            @test digits2(np) == vp
        end
    end
end