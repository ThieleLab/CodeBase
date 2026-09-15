function [A] = blkdiagi(B,j,n)
%% Description: 
%  return a block of size n and matrix B at the diagonal location j
%  BLKDIAG  Block diagonal concatenation of matrix input arguments.
%
%                                   |0 0 .. 0|
%   A = BLKDIAG(.,B,...)  produces  |0 B .. 0| of size n 
%                                   |0 0 ..  |
%
%   Class support for inputs:
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%      char, logical
%
%   See also DIAG, HORZCAT, VERTCAT
%
% USAGE:
%   [A] = blkdiagi(B,j,n)
%
% INPUT:
% B     Matrix
% n     size
%
% OUTPUT:
% A     block diagonal form
%
% .. Authors:
%       - Original author: Mohammad Faiz Khan.
%
% .. Last updated: 28/04/2023
%% Code
i=1;
A=[];
    while(i<n)
        if(i==j)
            A = blkdiag(A,B);
        elseif(i~=j)    
            A = blkdiag(A,0);
        end
    i=i+1;
    end
end