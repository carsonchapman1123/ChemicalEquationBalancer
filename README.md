# Chemical Equation Balancer

A Python script that can be used to balance chemical equations using regular expressions and SymPy. The balance function takes in a string representing a chemical equation with no leading coefficients.

For example,
```python
print(balance("CH4+O2=CO2+H2O"))
```
would output "CH4+2O2=CO2+2H2O" to the console.