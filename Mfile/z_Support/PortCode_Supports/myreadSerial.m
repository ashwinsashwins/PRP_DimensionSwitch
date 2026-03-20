function out = myreadSerial(s)
out=0;
if s.BytesAvailable
    out = fread(s,s.BytesAvailable);
end

