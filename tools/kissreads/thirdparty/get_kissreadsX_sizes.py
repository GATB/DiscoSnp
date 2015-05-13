#!/usr/bin/env python
import sys
f=open(sys.argv[1], "r")

while 1:
	com=f.readline() #>0__len__55707[0,791[
	if not com: break
	f.readline() # read the sequence.
	start=int(com.split('[')[1].split(',')[0].strip())
	stop=int(com.split('[')[1].split(',')[1].split('[')[0].strip())
	print stop-start
