export df_tlag, renamedf!
using DataFrames

function df_tlag(df::DataFrame, tlag::Int, spliter::String)
    df0 = copy(df)
    for i in 1:tlag
        df1 = copy(df)
        for j in 1:(i-1)
            delete!(df1, [size(df1)[1]])
        end
        renamedf!(df1, spliter * string(i))
        delete!(df0, [1])
        delete!(df1, [size(df1)[1]])

        df0 = hcat(df0, df1)
    end
    df0
end


function renamedf!(df::DataFrame, pt::String)
    nms = names(df)
    for i in 1:length(nms)
        nms[i] = nms[i] * pt
    end
    rename!(df, nms)
end