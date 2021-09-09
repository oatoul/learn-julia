using DataStructures

function gene_wise_consistency(X::Array{String}, Y::Array{Int})
    arr0 = String[]
    arr1 = String[]

    len = length(X)

    for i in 1:length(X)
        if Y[i] == 0
            push!(arr0, X[i])
        end

        if Y[i] == 1
            push!(arr1, X[i])
        end
    end

    cnt_arr0 = counter(arr0)
    cnt_arr1 = counter(arr1)
    f_cnt = 0

    for k in keys(cnt_arr0)
        f_cnt += min(cnt_arr0[k], cnt_arr1[k])
    end

    1 - f_cnt/len
end

# x = ["100", "001", "111", "100", "011"]
# y = [0, 1, 0, 1, 0]

# css = gene_wise_consistency(x, y)

# println(css)