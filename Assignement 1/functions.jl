# Functions for all recursion types
recursion64(a_n0::Float64, a_n1::Float64)       = Float64(2)*a_n0  - Float64(8)/Float64(9)  *a_n1;
recursion32(a_n0::Float32, a_n1::Float32)       = Float32(2)*a_n0  - Float32(8)/Float64(9)  *a_n1;
recursionBig(a_n0::BigFloat, a_n1::BigFloat)    = BigFloat(2)*a_n0 - BigFloat(8)/BigFloat(9)*a_n1;


# Plot the real value
function plot_true()
    x = 1:80
    y = 3/2*(2/3).^x
    plt[:plot](x,y, label="True Function", color="purple", ls="--", linewidth=3)
end
