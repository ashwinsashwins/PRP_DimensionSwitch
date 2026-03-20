function write_parallel_IO64(CONFIG,address,evcode)

%test for correct number of input arguments
if(nargin ~= 3)
    error('usage: write_parallel_IO64(CONFIG,address,data)');
end

io64(CONFIG.io.ioObj,address,evcode);
