function outputTable = getTable(traceTable, varName)
    C = cell(height(traceTable), 1);
    C(:) = {varName};
    rows = cellfun(@strcmp, traceTable.FileName, C);
    outputTable = traceTable(rows,:);
end