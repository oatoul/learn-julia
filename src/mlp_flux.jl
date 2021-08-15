export train_mlp_flux
using Flux
using Flux: gradient
using Flux.Optimise: update!
using DelimitedFiles, Statistics
using Parameters
using CSV
using DataFrames

# This replicates the housing data example from the Knet.jl readme. Although we
# could have reused more of Flux (see the mnist example), the library's
# abstractions are very lightweight and don't force you into any particular
# strategy.

# Struct to define hyperparameters
@with_kw mutable struct Hyperparams
    lr::Float64 = 0.1		# learning rate
    split_ratio::Float64 = 0.1	# Train Test split ratio, define percentage of data to be used as Test data
end

function get_processed_data(args, x_idx, y_idx)

    # println("X_idx $(x_idx)")
    # println("Y_idx $(y_idx)")

    df_origin = DataFrame(CSV.File("data/CDC15.tsv",drop=["Time"],type=Float32))
    df = Matrix(df_origin)'

    split_ratio = args.split_ratio # For the train test split

    # x = rawdata[1:13,:]
    # y = rawdata[14:14,:]
    x = df[x_idx,1:end-1]
    y = df[y_idx:y_idx,2:end]

    # Normalise the data
    x = (x .- mean(x, dims = 2)) ./ std(x, dims = 2)

    # Split into train and test sets
    split_index = floor(Int,size(x,2)*split_ratio)
    x_train = x[:,1:split_index]
    y_train = y[:,1:split_index]
    x_test = x[:,split_index+1:size(x,2)]
    y_test = y[:,split_index+1:size(x,2)]

    train_data = (x_train, y_train)
    test_data = (x_test, y_test)

    return train_data,test_data
end

# Struct to define model
mutable struct model
    W::AbstractArray
    b::AbstractVector
end

# Function to predict output from given parameters
predict(x, m) = m.W*x .+ m.b

# Define the mean squared error function to be used in the loss 
# function. An implementation is also available in the Flux package
# (https://fluxml.ai/Flux.jl/stable/models/losses/#Flux.Losses.mse).
meansquarederror(ŷ, y) = sum((ŷ .- y).^2)/size(y, 2)

function train_mlp_flux(x_idx::Array{Int64}, y_idx::Int64; kws...)

    # println("X_idx $(x_idx)")
    # println("Y_idx $(y_idx)")
    # Initialize the Hyperparamters
    args = Hyperparams(; kws...)
    
    # Load the data
    (x_train,y_train),(x_test,y_test) = get_processed_data(args, x_idx, y_idx)
    
    # The model
    m = model((randn(2,size(x_idx)[1])),[0.])
    
    loss(x, y) = meansquarederror(predict(x, m), y) 

    ## Training
    η = args.lr
    θ = params(m.W, m.b)

    for i = 1:500
        g = gradient(() -> loss(x_train, y_train), θ)
        for x in θ
            update!(x, g[x]*η)
        end
        if i%100==0
            @show loss(x_train, y_train)
        end
    end
    
    # Predict the RMSE on the test set
    err = meansquarederror(predict(x_test, m),y_test)

    err
end

# cd(@__DIR__)
# train([1,2,3], 4)