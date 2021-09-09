export dynamic_consistency
using DataStructures
using CSV
using DataFrames
import Statistics

function gene_wise_consistency(X::Vector, Y::Array{Int})
    arr0 = []
    arr1 = []

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


function get_conn_map(conns::Set)
    conn_map = Dict{String, Set{String}}()
    for con in conns
        r_t = split(con, "_")
        re = r_t[1]
        ta = r_t[2]
        vals = get!(Set{String}, conn_map, ta)
        push!(vals, re)
    end
    println(conn_map)
    conn_map
end


function dynamic_consistency(df::DataFrame, conns::Set)
    conn_map = get_conn_map(conns)
    gene_consis = []
    for (ta, res) in conn_map
        rea = toArray(res)
        X = copy(Tuple.(eachrow(df[!,rea])))
        Y = copy(df[!,ta])
        # t + 1
        pop!(X)
        popfirst!(Y)

        g_c = gene_wise_consistency(X, Y)
        println("$(ta) $(g_c)")
        push!(gene_consis, g_c)
    end
    Statistics.mean(gene_consis)
end

function toArray(s::Set)
    (x -> x).(s)
end

# x = ["100", "001", "111", "100", "011"]
# y = [0, 1, 0, 1, 0]

# css = gene_wise_consistency(x, y)

# println(css)

function get_expect_CDC15()
    res = Set()
    push!(res, "FKH2_SWI5")
    push!(res, "FKH2_ACE2")

    push!(res, "ACE2_CLN3")

    push!(res, "SWI6_NDD1")
    push!(res, "SWI6_SWI4")

    push!(res, "NDD1_ACE2")
    push!(res, "NDD1_SWI5")

    push!(res, "SWI5_CLN3")

    push!(res, "MCM1_ACE2")
    push!(res, "MCM1_SWI5")
    push!(res, "MCM1_CLN3")
    push!(res, "MCM1_SWI4")

    push!(res, "CLN3_SWI4")
    push!(res, "CLN3_MBP1")

    push!(res, "MBP1_NDD1")
    push!(res, "MBP1_SWI4")

    push!(res, "SWI4_NDD1")

    res
end

function get_NNBNI_CDC15()
    res = Set()
    push!(res, "FKH2_SWI5")
    push!(res, "FKH2_ACE2")

    push!(res, "ACE2_CLN3")
    push!(res, "ACE2_MCM1")
    push!(res, "ACE2_SWI4")
    push!(res, "ACE2_SWI5")
    push!(res, "ACE2_SWI6")

    push!(res, "SWI6_NDD1")
    push!(res, "SWI6_SWI4")

    push!(res, "NDD1_ACE2")
    push!(res, "NDD1_SWI5")
    push!(res, "NDD1_CLN3")
    push!(res, "NDD1_MBP1")

    push!(res, "SWI5_CLN3")

    push!(res, "MCM1_ACE2")
    push!(res, "MCM1_SWI5")
    push!(res, "MCM1_NDD1")
    push!(res, "MCM1_CLN3")
    push!(res, "MCM1_SWI4")

    push!(res, "CLN3_SWI4")
    push!(res, "CLN3_MBP1")
    push!(res, "CLN3_MCM1")

    push!(res, "MBP1_FKH2")
    push!(res, "MBP1_SWI6")
    push!(res, "MBP1_NDD1")
    push!(res, "MBP1_SWI4")

    push!(res, "SWI4_ACE2")
    push!(res, "SWI4_NDD1")

    res
end

# df_origin_bool = DataFrame(CSV.File("data/CDC15_bool.tsv",drop=["Time"],type=Int))
# println(df_origin_bool)

# dc = dynamic_consistency(df_origin_bool, get_NNBNI_CDC15())
