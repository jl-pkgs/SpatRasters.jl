using SpatRasters, ArchGDAL, NetCDFTools
using Ipaper, DataFrames, RTableTools, JLD2, NaNStatistics

# # xlim = c(109.4, 111.6), ylim = c(31.2, 33.4)
# b = bbox(109.4, 31.4, 111.6, 33.4)
b = bbox(109.4, 31.2, 111.6, 33.4)
target = make_rast(; b, cellsize=1 / 120 * 2)

st = fread("./Project_十堰/data/十堰_雨量站_sp237_v20250824.csv")
df = fread("./Project_十堰/data/十堰_prcp_sp237_mat_v20250824.csv")

# 5.45%数据缺失
dates = DateTime.(df.time, "yyyy-mm-ddTHH:MM:SSZ")
data = Matrix(df[:, 2:end])' |> drop_missing |> collect

# 无降水不插值
_prcp = NaNStatistics.nansum(data, dims=1)[:]
inds = findall(_prcp .>= 1.0) # 所有站点总降水大于1mm日期，才进行插值
info = DataFrame(; time=dates[inds], inds, prcp=_prcp[inds])

## 部分数据存在空值，如何处理
X = [st.lon st.lat] # d.alt
Y = data[:, inds]
points = map(x -> (x[1], x[2]), eachrow(X))

nlon, nlat = size(target)[1:2]
ntime = length(dates)


# TPS: 不能存在NaN，因此采用IDW插值的结果进行填补
if false
  locs = st_location(target, points)
  _inds = map(p -> LinearIndices(target.A)[p...], locs) # locs2映射到2维

  Z = nc_read("ShiYan_Prcp_Gauged237_201404-202501_2km_IDW(r=30,nmax=6).nc")
  _Z = reshape(Z, nlon * nlat, :)[_inds, :]

  data2 = deepcopy(data)
  inds_bad = isnan.(data)
  data2[inds_bad] .= _Z[inds_bad]

  jldsave("ShiYan_Pobs_interpolated_by_IDW.jld2", true; st, dates, P=data2)
end

# 缺失比例降为0.25%，余下的部分赋值为0
data2 = jldopen("./Project_十堰/data/ShiYan_Pobs_interpolated_by_IDW.jld2")["P"]
Y2 = data2[:, inds]
Y2[isnan.(Y2)] .= 0.0 # IDW插值，有些区域无法覆盖


# TPS
for λ in [0.001, 0.01]
  @time ra_tps1 = interp_tps(X, Y, target; λ)

  R = zeros(Float32, size(target)[1:2]..., length(dates))
  R[:, :, inds] .= ra_tps1.A

  lon, lat = st_dims(target)
  dims = (; lon, lat, time=dates)
  fout = "ShiYan_Prcp_Gauged237_201404-202501_2km_TPS(λ=$λ).nc"
  ncsave(fout, true, (; units="mm h-1"); dims, P=R)
end


# IDW
if false
  radius = 30
  nmax = 6
  @time ra_idw = interp(X, Y, target; method="idw", radius, nmax, m=2)

  R = zeros(Float32, size(target)[1:2]..., length(dates))
  R[:, :, inds] .= ra_idw.A

  lon, lat = st_dims(target)
  dims = (; lon, lat, time=dates)
  fout = "ShiYan_Prcp_Gauged237_201404-202501_2km_IDW(r=30,nmax=6).nc"
  ncsave(fout, true, (; units="mm h-1"); dims, P=R)
end

# ADW
if false
  radius = 30
  nmax = 6

  for m in [2, 4, 6]
    ra = interp(X, Y, target; method="adw", radius, cdd=25, nmax, do_angle=true, m)

    R = zeros(Float32, size(target)[1:2]..., length(dates))
    R[:, :, inds] .= ra.A

    lon, lat = st_dims(target)
    dims = (; lon, lat, time=dates)
    fout = "ShiYan_Prcp_Gauged237_201404-202501_2km_ADW(r=30,nmax=6,m=$m,cdd=25).nc"
    ncsave(fout, true, (; units="mm h-1"); dims, P=R)
  end
end
