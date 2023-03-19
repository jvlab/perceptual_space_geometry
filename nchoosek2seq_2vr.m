function z=nchoosek2seq_2vr(a,amax)
% z=nchoosek2seq_2vr(a,amax) maps an increasing sequence of 2 natural numbers a into an integer, using lexicographical ordering 
% with rightmost value of a increasing first
%
%  z=1,2,3 corresponds to a=[1 2],[1 3],[1 4],..., [1 amax],[2 3], [2 4]
% contrast with nchoosek2seq_2v, for which z=1,2,3 corresonds to a=[1 2],[1 3],[2 3],[2 4]; ...
%
% a: a sequence of increasing triplets of integers, or a row of such sequences
% amax: maximum possible value of a
% 
%   See also:  SEQ2NCHOOSEK, NCHOOSEK2SEQ_2V, NCHOOSEK2SEQ_VR.
acombs=nchoosek(amax,2);
z=acombs+1-nchoosek2seq_2v(amax+1-fliplr(a));
return
