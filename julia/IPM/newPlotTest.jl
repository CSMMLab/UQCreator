using PyCall
@pyimport matplotlib.animation as anim
using PyPlot

function showgif(filename)
    open(filename) do f
        base64_video = base64encode(f)
        display("text/html", """<img src="data:image/gif;base64,$base64_video">""")
    end
end

fig = figure()
ax = axes()

x = [-1.5:0.05:1.0;]
y = [-1.0:0.05:1.0;]
z = x' .+ im*y
w0 = sqrt.(z.+1)./z

# plotw has period 1.
# plotw(t+1) == plotw(1)
function plotw(t)
    clf()
    title("\$t = $t\$")
    xlim(-1.4,0.9)
    ylim(-0.9,0.9)
    w = exp(pi*im*t)*w0
    streamplot(x,y,real(w),-imag(w))
    axes()[:set_aspect]("equal")
end

n=10 # frame/loop
l=1 # loop
interval=200 # milli seconds

withfig(fig) do
    for k in 1:n*l
        plotw(k/n)
        savefig(@sprintf("test5_%04d",k), bbox_inches="tight")
    end
end

#run(`convert -delay $(interval/10) -loop 0 tmp0\*.png tmp.gif`)
run(`magick -delay $(interval/10) -loop 0 test5_\*.png test5.gif`)

showgif("test5.gif")