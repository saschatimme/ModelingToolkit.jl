using ModelingToolkit

l4(x,y,z,w,a,b,c,d) = x*a+y*b+z*c+w*d
register_derivatives(l2, 8)
@register dump()
dump()
dot_2(x,a) = x*a

@register_derivatives

expand_derivatives(Differential(x[1])(dot_2(x...,a...)))
@macroexpand @register dot_1(x,y,z,w,a,b,c,d)

@macroexpand (ModelingToolkit.@register dot_2(x,d))

@variables x[1:4] a[1:4]


expand_derivatives(Differential(x[1])(dot_1(x...,a...)^2))


function register_with_derivatives(expr)
    name = expr.args[1]
    n = length(expr.args) - 1
    op = @eval begin
        let
            ModelingToolkit.@variables x[1:$n]
            $name(x...)
        end
    end
    ModelingToolkit.@variables x[1:n]
    var_names = simplified_expr.(x)
    esc(Expr(:block, map(1:n) do i
        quote
            function ModelingToolkit.derivative(::typeof($name), args::NTuple{$n,Any}, ::Val{$i})
                $(Expr(:tuple, var_names...)) = args
                $(simplified_expr(expand_derivatives(Differential(x[i])(op))))
            end
        end
    end..., :(@register $expr)))
end

macro register_with_derivatives(expr)
    register_with_derivatives(expr)
end
register_with_derivatives(:())

@macroexpand(@register_with_derivatives l4(x,y,z,w,a,b,c,d))

@register_with_derivatives l4(x,y,z,w,a,b,c,d)

ModelingToolkit.@variables x[1:4] a[1:4]

g = 1 / l4(x..., a...)

e1 = Differential(x[1])(g)
expand_derivatives(e1)
e = Differential(x[1])()
expand_derivatives(e)

x₁, x₂, x₃, x₄ = x
a₁, a₂, a₃, a₄ = a

f(z) = z
ModelingToolkit.derivative(::typeof(f), ::Tuple{Any}, ::Val{1}) = 1
@register f(z)
@variables z w
expr = ModelingToolkit._simplify_constants((-1) * (1 / z) * z, true)
expand_derivatives(Differential(z)(expr))

expand_derivatives(Differential(z)(expr))
expr = (-1 * ((1 / l4(x₁, x₂, x₃, x₄, a₁, a₂, a₃, a₄)) / l4(x₁, x₂, x₃, x₄, a₁, a₂, a₃, a₄))) * a₁
expr.args
ModelingToolkit._simplify_constants(expr, true)

ModelingToolkit.hessian(g, x)

op = dot_2(x...)
Differential()

register_derivatives(dot_2, 2)

@register dot_2(x,a)


@variables x a


show(simplified_expr(expand_derivatives(Differential(x)(1/dot_2(x,a))) ))

f(x,a) = -1 * ((1 / dot_2(x, a)) / dot_2(x, a)) * a

@code_native f(2.3,4.2)
function _register_derivaties(E::Expr)
    name = E.args[1].args[1]
    arguments = E.args[1].args[2:end]
    N = length(arguments)
    vars = [Variable(arguments[i])() for i in 1:N]
    symbolic_output = @eval begin
        let
            $(arguments...) =
        end
    end
    map(1:N) do i
        quote
            function ModelingToolkit.derivative(typeof($name), args::NTuple{$N,Any}, ::Val{$i})
                $()
            end
        end
    end
end


E = :(function dot_1(x,y,z,w,a,b,c,d)
        x*a+y*b+z*c+d*w
end)

_register_derivaties(E)
E.args[1].args

@variables args[1:8]
vars = map(x -> Variable(x)(), E.args[1].args[2:end])
x, y, z, w, a, b, c, d = vars
body = eval(E.args[2])



map(enumerate(vars)) do ((i,var),)
    @eval function ModelingToolkit.derivative(::typeof($(E.args[1].args[1])), args::NTuple{$(length(vars)),Any}, ::Val{$i})
            x, y, z, w, a, b, c, d = args
            simplified_expr(expand_derivatives(Differential(x)(body)))
        end
    end
end
simplified_expr(expand_derivatives(Differential(x)(body)))
