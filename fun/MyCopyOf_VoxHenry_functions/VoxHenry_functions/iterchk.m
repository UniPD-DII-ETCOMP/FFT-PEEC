function [atype,afun,afcnstr] = iterchk(A)
%ITERCHK  Checks arguments to iterative methods.
%   [ATYPE,AFUN,AFCNSTR] = ITERCHK(A) returns the following:
%   ATYPE is either 'matrix', 'function', 'expression' or 'inline object'.
%   AFUN is the function name or inline object.
%   AFUN is '' if ATYPE is 'matrix'.
%   AFCNSTR is the function name if ATYPE is 'function'.
%   AFCNSTR is the formula of the function if ATYPE is 'expression' or
%   'inline object'.  AFCNSTR is '' if ATYPE is 'matrix'.
%
%   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ.

%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 1.8.4.2 $ $Date: 2004/12/06 16:35:56 $

% if running under Octave, use the fcnchk_octave() alternative version
if(exist ('OCTAVE_VERSION', 'builtin') > 0)
   [afun,afunmsg] = fcnchk_octave(A); 
else
   [afun,afunmsg] = fcnchk(A);
end

% 'afunmsg' contains any error message from fcnchk
if isempty(afunmsg)
   if isa(afun,'inline')      
      if isa(A,'inline')
         atype = 'inline object';
      else
         atype = 'expression';
      end
      afcnstr = formula(afun);
   else % both function_handles @fun and function names 'fun'
      atype = 'function';
      if isa(A,'function_handle')
          afcnstr = func2str(A);
      else
          afcnstr = A;
      end
   end
elseif isa(A,'float')
   afun = A;
   atype = 'matrix';
   afcnstr = '';
else
   error('MATLAB:iterchk:InvalidInput',...
         'Argument must be a floating point matrix or a function handle.');
end
