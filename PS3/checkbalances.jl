using DelimitedFiles

A=readdlm("atom.dat")
display(A)

S=readdlm("Network.dat")
display(S)

e=transpose(A)*S
display(e)
