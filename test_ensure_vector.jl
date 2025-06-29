# 测试 ensure_vector 函数
include("src/SpatRaster.jl")

# 测试 ensure_vector 函数
println("测试 ensure_vector 函数:")

# 测试 nothing
result = ensure_vector(nothing)
println("ensure_vector(nothing) = ", result)

# 测试标量
result = ensure_vector("band1")
println("ensure_vector(\"band1\") = ", result)

# 测试已经是 vector
result = ensure_vector(["band1", "band2"])
println("ensure_vector([\"band1\", \"band2\"]) = ", result)

# 测试数字标量
result = ensure_vector(42)
println("ensure_vector(42) = ", result)

# 测试 ensure_bands_time_vector 函数
println("\n测试 ensure_bands_time_vector 函数:")

bands, time = ensure_bands_time_vector("band1", 123)
println("ensure_bands_time_vector(\"band1\", 123) = ", (bands, time))

bands, time = ensure_bands_time_vector(["band1", "band2"], nothing)
println("ensure_bands_time_vector([\"band1\", \"band2\"], nothing) = ", (bands, time))

bands, time = ensure_bands_time_vector(nothing, [1, 2, 3])
println("ensure_bands_time_vector(nothing, [1, 2, 3]) = ", (bands, time))
