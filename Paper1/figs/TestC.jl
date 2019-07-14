using PyPlot

NMax = 50;

entries = zeros(NMax);
for N = 1:NMax
    #C = zeros(N,N,N);
    counterNonZero = 0;
    for i = 1:N
        for j = 1:N
            for k = 1:N
                s = 0.5*(i+j+k)
                if s >= i && s >= j && s >= k
                    counterNonZero += 1;
                    #C[i,j,k] = factorial(i)*factorial(j)*factorial(k)/(factorial(s-i)*factorial(s-j)*factorial(s-k));
                else
                    #C[i,j,k] = 0;
                end
            end
        end
    end
    entries[N] = counterNonZero;
end

fig, ax = subplots(figsize=(15, 8), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax.plot((1:NMax).^3,entries, "r--o", linewidth=2, alpha=0.6)
ax.set_xlabel(L"$N^3$", fontsize=20);
ax.set_ylabel("entries in C", fontsize=20);
ax.tick_params("both",labelsize=20)
ax.set_xlim([0,NMax^3])
ax.set_ylim([0,entries[end]+100])