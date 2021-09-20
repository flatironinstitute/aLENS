# Converts MT number density (1/um^3)
# to tubulin dimer (1 dimer = 1 alpha + 1 beta)
# concentration umol/L
# 1625 Tubulin dimer / um

import argparse

parser = argparse.ArgumentParser(
    description='Convert MT number density to tubulin dimer molar concentration.')
parser.add_argument('N', type=int, help='MT number')
parser.add_argument('L', type=float, help='MT average length')
parser.add_argument('-d', type=float, dest='d',
                    help='MT diameter', default=0.025)
parser.add_argument('x', type=float, help='box x in um')
parser.add_argument('y', type=float, help='box y in um')
parser.add_argument('z', type=float, help='box z in um')

args = parser.parse_args()

number = args.N
length = args.L
x = args.x
y = args.y
z = args.z
radius = args.d/2

print(args)

volume = x*y*z
volOneMT = 3.1416*radius**2*length

print("tubulin dimer concentration in umol/L: ", 1625 *
      length*number/(volume)*(10**15)/(6.02*10**23)*10**6)

print("MT volume fraction: ", volOneMT*number/volume)
