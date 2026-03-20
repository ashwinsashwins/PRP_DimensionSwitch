function [port_state]=read_parallel(portaddress);
%   Usage of read_parallel(portaddress)
%   copy inpout32.dll to c:\windows\system32\
%   include parallel_in.dll in matlab path
%   usage:
% read_parallel(888): portaddress-decimal(888) to read from
% read_parallel('378'): portaddress-hexadecimal(x378) to read from
% read_parallel : default portaddress-decimal(888) to read from

if (nargin == 0)
portaddress=888;
elseif (nargin == 1)
if ischar(portaddress)
  portaddress=hex2dec(portaddress);  
end
end
 port_state=parallel_in(portaddress);
 