#my first try at input my own function into Python

def fib(arg):
	"Prints a Fibonacci series up to the value of the input arg"
	a, b = 0, 1
	while b < arg:
		print b,
		a, b = b, a + b
