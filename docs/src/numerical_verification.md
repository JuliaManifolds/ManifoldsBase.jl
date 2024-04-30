# Numerical Verification

```@autodocs
Modules = [ManifoldsBase]
Pages = ["src/numerical_checks.jl"]
Order = [:macro, :type, :function]
Private = false
Public = true
```

## Internal functions

The following functions split the check into several parts, for
example looking for the best fitting window and finding out the best slope,
or plotting the slope.

```@autodocs
Modules = [ManifoldsBase]
Pages = ["src/numerical_checks.jl"]
Order = [:macro, :type, :function]
Private = true
Public = false
```