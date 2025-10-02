#python example_doublepowerlaw.py density="potential_spherical.ini" > spherical.log
#python example_doublepowerlaw.py density="potential_oblate_0.5.ini" > oblate_0.5.log
python example_doublepowerlaw.py density="potential_oblate_0.5.ini" nbody=100000 out=1e5_oblate_0.5.dat => oblate_0.5_1e5.log

