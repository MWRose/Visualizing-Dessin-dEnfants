using HomotopyContinuation
using Plots
using DelimitedFiles
using DataFrames
using CSV
using SymPy
using JSON
using HTTP


@var x
h = 0.0001

function get_roots(f)
    @var y
    roots = []
    # define the polynomials
    F = System([f])
    result = HomotopyContinuation.solve(F)
    for root in results(result)
        push!(roots, root.solution[1])
    end
    return roots
end

function get_data(name)
    url = string("https://beta.lmfdb.org/Belyi/download_galmap_to_text/", name)
    r = String(HTTP.request("GET", url).body)
    brace_1 = findfirst("{", r)
    brace_2 = findfirst("}", r)
    data_string = r[brace_1[1]:brace_2[1]]
    data = JSON.parse(data_string)
    return data
end

function get_functions(data)

    base_map = data["map"]
    b_string = replace(base_map, "nu" => data["embeddings"][1][1])
    # b_string = base_map
    b_sym = simplify(sympify(b_string))
    println(b_sym)

    n = eval(convert(Expr, numer(b_sym)))
    d = eval(convert(Expr, denom(b_sym)))
    belyi = b_sym

    body = convert(Expr, belyi)
    syms = Symbol.(free_symbols(belyi))
    bt = eval(Expr(:function, Expr(:call, gensym(), syms...), body))

    ex = diff(belyi)
    bodydiff = convert(Expr, ex)
    symsdiff = Symbol.(free_symbols(ex))
    dt = eval(Expr(:function, Expr(:call, gensym(), syms...), bodydiff))

    return belyi, n, d, bt, dt

end

function get_functions_from_string(b)
    b_string = b
    # b_string = base_map
    b_sym = simplify(sympify(b_string))
    println(b_sym)

    n = eval(convert(Expr, numer(b_sym)))
    d = eval(convert(Expr, denom(b_sym)))
    belyi = b_sym

    body = convert(Expr, belyi)
    syms = Symbol.(free_symbols(belyi))
    bt = eval(Expr(:function, Expr(:call, gensym(), syms...), body))

    ex = diff(belyi)
    bodydiff = convert(Expr, ex)
    symsdiff = Symbol.(free_symbols(ex))
    dt = eval(Expr(:function, Expr(:call, gensym(), syms...), bodydiff))

    return belyi, n, d, bt, dt
end

function rk4(start, stop, initialQ, initialE, dt)
    epoints = []
    qs = []
    rs = []
    q = initialQ
    e = initialE
    push!(epoints, e)

    while q < stop
        nexth = 1/2 * abs(dt(e)) *(1 + abs(e)^2) * h 
        nextq = q + nexth
        
        
        k1 = nexth * (1/dt(e))
        k2 = nexth * (1/dt(e + k1/2))
        k3 = nexth * (1/dt(e + k2/2))
        k4 = nexth * (1/dt(e + k3))
        
        nexte = e + (1/6) * (k1 + 2*k2 + 2*k3 + k4)

        if isnan(nexte)
            println("NEXT E NAN: ", e)
        end

        push!(epoints, nexte)
        q = nextq
        e = nexte
    end
    q = initialQ
    e = initialE
    
    while q > start
        nexth = 1/2 * abs(dt(e)) *(1 + abs(e)^2) * h 
        nextq = q - nexth

        k1 = nexth * (1/dt(e))
        k2 = nexth * (1/dt(e - k1/2))
        k3 = nexth * (1/dt(e - k2/2))
        k4 = nexth * (1/dt(e - k3))
        
        nexte = e - (1/6) * (k1 + 2*k2 + 2*k3 + k4)

        if isnan(nexte)
            println("NEXT E NAN: ", e)
        end
        
        push!(epoints, nexte)
        q = nextq
        e = nexte
        
    end
    return epoints, qs, rs
end

function stereographic(c)
    x = (2 * real(c)) / (abs(c) ^ 2 + 1)
    y = (2 * imag(c)) / (abs(c) ^ 2 + 1)
    z = (abs(c)^2 - 1) / (abs(c)^2 + 1)
    if isnan(x)
        println("THE VALUE WAS NAN: ", c)
    end
    return (x, y, z)
end

function compute_stereo(points)
    xs = []
    ys = []
    zs = []
    for point in points
        x, y, z = stereographic(point)
        x = round(Float16(x, RoundDown), digits=3)
        y = round(Float16(y, RoundDown), digits=3)
        z = round(Float16(z, RoundDown), digits=3)

        push!(xs, x)
        push!(ys, y)
        push!(zs, z)
    end
    return xs, ys, zs
end

function compute_colors(points, bt)
    colors = []
    for point in points
        color = real(bt(point))
        if color > 1
            color = 1
        end
        if color < 0
            color = 0
        end
        r = 12
        color = (1 / (2 * r)) * (log(2.75, (1 + color * (exp(r) - 1))/(1 + color * (exp(-1 * r) - 1))))
        color = round(Float16(color, RoundDown), digits=3)
        push!(colors, color)
    end
    return colors
end

function create_pc(xs, ys, zs, colors, sampling, name)
    sampling = 1
    l = length(xs)/sampling
    io = open(name, "w"); 
    # write(io,
    # "
    # # .PCD v0.7 - Point Cloud Data file format
    # VERSION 0.7
    # FIELDS x y z
    # SIZE 4 4 4
    # TYPE F F F
    # COUNT 1 1 1
    # WIDTH $(length(xs)/sampling)
    # HEIGHT 1
    # VIEWPOINT 0 0 0 1 0 0 0
    # POINTS $(length(xs)/sampling)
    # DATA ascii\n
    # ")
        for i in 1:(length(xs))
            write(io, "$(xs[i]) $(ys[i]) $(zs[i]) $(colors[i]) $(colors[i]) $(colors[i])\n")
        end
    
    close(io);
end

function main()
    # name = "8T49-5.3_5.3_3.1.1.1.1.1-a"
    # name = "8T43-6.1.1_6.1.1_2.2.2.2-a"
    name = "7T7-4.2.1_4.3_2.2.2.1-a"
    # name = "((x^(20) + 228*x^(15) + 494*x^(10) - 228*x^5 + 1)^3 ) / (1728 * x^5 * (x^(10) - 11*x^5 - 1)^5 )"
    # name = "(1728 * x^5 * (x^(10) - 11*x^5 - 1)^5 ) / (x^(20) + 228*x^(15) + 494*x^(10) - 228*x^5 + 1)^3"
    # name = "9T33-5.3.1_2.2.2.2.1_3.3.3-a"
    # name = "x^5"
    outfile = "newPC.txt"
    data = get_data(name)
    
    belyi, n, d, bt, dt = get_functions(data)
    # belyi, n, d, bt, dt = get_functions_from_string(name)
    println("Finding Roots")
    roots = get_roots(n - 1/2 * d)
    
    println("Finished\nROOTS FOUND: $(length(roots))")

    println("Computing edges")
    rs = []
    qs = []
    points = []
    for root in roots
        edge = Base.invokelatest(rk4, 0, 1, 0.5, root, dt)
        points = vcat(points, edge[1])
        qs = vcat(qs, edge[2])
        rs =vcat(rs, edge[3])
    end
    println("Finished\nPOINTS FOUND: $(length(points))")

    println("Projecting to sphere")
    xs, ys, zs = compute_stereo(points)
    println("Finished")
    println("Calculating Colors")
    colors = Base.invokelatest(compute_colors, points, bt)
    println("Finished")
    println("Writting data to file")
    create_pc(xs, ys, zs, colors, 1, outfile)
    println("Finsihed")
    println("FILE: $(outfile)")
end

main()

