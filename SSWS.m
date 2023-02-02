function E = SSWS(c, F, k)

c = cell2mat(c);

c = reshape(c.', 1, 1, 1, 3, []);

E = -1i.*k.*sum(c.*cat(5, F{:}), 5);

end