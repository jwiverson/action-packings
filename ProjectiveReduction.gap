ProjectiveReduction:=function(G)
	#Input:		A Gram matrix G of some frame
	#Output:	The Gram matrix of a projective reduction of the same frame

	local n,unique,i,row,fn,new,j,oldRow,phase;

	n:=Size(G);
	unique:=[];

	for i in [1..n] do
		#check to see if the i-th row is a multiple of one we've seen
		row:=G[i];;
		fn:=PositionProperty(row,x->x<>0); #position of first nonzero entry
	
		new:=true; #it will be false if it matches an old row
		for j in unique do
			oldRow:=G[j];
			if oldRow[fn] = 0 then
				#they definitely don't match
				continue;
			else
				phase:=row[fn]/oldRow[fn];
				if row = phase*oldRow then
					#we've seen this one
					new:=false;
					break;
				fi;
			fi;
		od;

		if new then
			AddSet(unique,i);
		fi;
	od;

	return G{unique}{unique};
end;