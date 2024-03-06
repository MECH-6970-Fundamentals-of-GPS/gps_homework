function [CA] = genCA(sv, bits)
    tap = [2 6;
           3 7;
           4 8;
           5 9;
           1 9;
           2 10;
           1 8;
           2 9;
           3 10;
           2 3;
           3 4;
           5 6;
           6 7;
           7 8;
           8 9;
           9 10;
           1 4;
           2 5;
           3 6;
           4 7;
           5 8;
           6 9;
           1 3;
           4 6;
           5 7;
           6 8;
           7 9;
           8 10;
           1 6;
           2 7;
           3 8;
           4 9;
           5 10;
           4 10;
           1 7;
           2 8;
           4 10];

    if ~exist('bits','var')
        bits = 1023;
    end

    G1 = ones(1,10);
    G2 = ones(1,10);
    CA = zeros(1,bits);
    for i = 1:bits
        newBit1 = bitxor(G1(10),G1(3));
        newBit2 = bitxor(G2(2), bitxor(G2(3), bitxor(G2(6), bitxor(G2(8), bitxor(G2(9),G2(10))))));
        G2i = bitxor(G2(tap(sv,1)), G2(tap(sv,2)));
        CA(i) = bitxor(G1(10),G2i);
        G1 = [newBit1, G1(1:9)];
        G2 = [newBit2, G2(1:9)];
    end
end
