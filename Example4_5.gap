#This file contains the code behind Example 4.5
#of the paper "Optimal Line Packings From Finite Group Actions",
#by Joseph W. Iverson, John Jasper, and Dustin G. Mixon.
#
#The file "SphericalFunctions.gap" should be available through 
#the same source as this file.


#Create a method to search through a table of spherical functions
#and return the sizes of any ETFs that can be made using projective
#reduction.

search:=function(tbl)
	#Input:		A table of spherical functions, as returned by
	#		SphericalFunctions(G,H)
	#Output:	This function returns a list [complex,real],
	#		where complex is a list of all sizes of complex
	#		ETFs that can be made from rows of this table
	#		after projective reduction, and real is a similar
	#		list for real ETFs.
	#Note:		If there are more than 20 rows to check, we don't attempt it,
	#		and merely return [].

	local fullGroupSize,d,newTbl,scaledTbl,i,j,row,keepers,it,v,sum,rank,absSqr,stabPos,stabSize,n,real,complex,r,z;

	fullGroupSize:=Sum(tbl[1]);

	d:=Length(tbl);

	#Cut out anything with multiplicity greater than one
	newTbl:=[tbl[1]];
	for i in [2..d] do
		if tbl[i][1]=1 then
			Add(newTbl,tbl[i]);
		fi;
	od;
	tbl:=newTbl;
	d:=Length(tbl);

	if d>21 then
		#there are too many possibilities to reasonably check!
		return [];
	fi;

	#Scale the bottom rows by their multiplicities
	scaledTbl:=[tbl[1]];
	for i in [2..d] do
		row:=[];
		for j in [1..Size(tbl[i])-1] do
			Add(row,tbl[i][j]*tbl[i][Size(tbl[i])]);
		od;
		Add(row,tbl[i][Size(tbl[i])]);
		Add(scaledTbl,row);
	od;

	it:=IteratorOfCombinations([2..d]);

	real:=[];
	complex:=[];

	for v in it do
		if Length(v) = 0 then
			#Ignore the empty set
			continue;
		fi;

		sum:=Sum(v,i->scaledTbl[i]); #the sum of the rows indexed by v
		rank:=sum[Size(sum)];
		if rank = 1 then
			#We don't want any 1-dimensionals
			continue;
		fi;

		absSqr:=List([1..Size(tbl[i])-1],j->sum[j]*ComplexConjugate(sum[j])); 
		#the squares of the absolute values of the angles

		#Make sure it's not a simplex or an ONB
		stabPos:=Positions(absSqr,rank^2); #the positions corresponding to the projective stabilizer
		stabSize:=Sum(stabPos,j->tbl[1][j]); #the size of the projective stabilizer
		n:=fullGroupSize/stabSize; #the number of vectors in the frame corresponding to v
		if n=rank or n=rank+1 then
			#At best, we get an ONB or a simplex
			continue;
		fi;

		if Size(AsSet(absSqr))=2 then
			#We found a good one!

			#Check to see if it's real
			r:=true;
			for z in sum do
				if ComplexConjugate(z)<>z then
					#One of the inner products is not real
					r:=false;
					break;
				fi;
			od;

			#Record the result
			if r then
				#It's a real ETF
				AddSet(real,[rank,n]);
			else
				#It's a complex ETF
				AddSet(complex,[rank,n]);
			fi;
		fi;
	od;

	return [complex,real];
end;


#Create a list of transitive permutation groups
#of degree no more than 30 and not equal to 28, and of size no more
#than 100,000.

all:=AllTransitiveGroups(NrMovedPoints,Concatenation([1..27],[29,30]),Size,[1..100000]);;


#Make dictionaries to record the groups that produced ETFs of various sizes, real and complex,
#and lists to record the keys of those dictionaries

rDict:=NewDictionary([1,1],true); #keys are sizes of real ETFs
cDict:=NewDictionary([1,1],true); #keys are sizes of complex ETFs
rSizes:=[];; #list of keys to rDict
cSizes:=[];; #list of keys to cDict


#Iterate through the permutation groups, checking transitivity levels 1 and 2, and
#recording all ETF sizes made with the dictionaries.
#
#WARNING: This takes a very long time, on the order of several days.

Read("SphericalFunctions.gap");

for G in all do
	T:=Minimum(2,Transitivity(G));
	for t in [1..T] do
		#Create the stabilizer H for the action of G on t points
		x0:=[];
		for i in [1..t] do
			Add(x0,i);
		od;

		H:=Stabilizer(G,x0,OnTuples);

		#Compute the table of spherical functions and look for ETFs
		tbl:=SphericalFunctions(G,H);
		s:=search(tbl);
		
		#Get the group info if it made anything
		if Size(s[1])>0 or Size(s[2])>0 then
			m:=NrMovedPoints(G);
			k:=TransitiveIdentification(G);
		fi;

		#Record the results in the dictionaries
		for size in s[1] do
			#Use the complex dictionary
			if size in cSizes then
				#We have already seen this size
				grps:=LookupDictionary(cDict,size);
				AddSet(grps,[m,t,k]);
				AddDictionary(cDict,size,grps);
			else
				#This is a new size
				AddSet(cSizes,size);
				AddDictionary(cDict,size,[[m,t,k]]);
			fi;
		od;

		for size in s[2] do
			#Use the real dictionary
			if size in rSizes then
				#We have already seen this size
				grps:=LookupDictionary(rDict,size);
				AddSet(grps,[m,t,k]);
				AddDictionary(rDict,size,grps);
			else
				#This is a new size
				AddSet(rSizes,size);
				AddDictionary(rDict,size,[[m,t,k]]);
			fi;
		od;
	od;
od;


#Print the results
List(rSizes,size->List(LookupDictionary(rDict,size),info->Concatenation(size,info)));
List(cSizes,size->List(LookupDictionary(cDict,size),info->Concatenation(size,info)));

