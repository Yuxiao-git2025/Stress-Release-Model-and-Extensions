function params=LSRM_PackingParam(A,B,C)
    params=[A;B;reshape(C,[],1)];       % row-major packing
    % Notice the dimension is adjusted to fit the unpack procedure
end