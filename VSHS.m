function F = VSHS(e, m, E, M)

F = sum(reshape(e, 1, 1, 1, []).*cat(4, E{:}) + reshape(m, 1, 1, 1, []).*cat(4, M{:}), 4);

end