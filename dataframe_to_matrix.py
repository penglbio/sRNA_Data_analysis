#!/usr/bin/env python
import sys

k={}
win_list=[]
len_list=[]
with open(sys.argv[1]) as f:
	for line in f:
		(win_id,length,hna)=line.strip().split("\t")
		if win_id in win_list:
			len_list.append(length)
			k[length]=line.strip()
		else:
			win_list.append(win_id)
			if len(len_list) != 0:
				for i in range(18,61):
					i=str(i)
					if i in len_list:
						print(k[i])
					else:
						print('%s\t%s\t%s' % (win_list[-2],i,0.00))
				len_list=[]
				k={}
			else:
				len_list.append(length)	
				k[length]=line.strip()
