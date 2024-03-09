function saveModel(result::SuccessResult, file_name::String; overwrite = true)
  data = result.data
  sol = result.solution
  data_dfs = data2DataFrame(data)
  sol_dfs = solution2DataFrame(sol)
  sol_sheet_names = keys(sol_dfs)
  XLSX.writetable(
    "$file_name.xlsx",
    ["solution_$(String(name))" => sol_dfs[name] for name in sol_sheet_names]...,
    [String(name) => data_dfs[name] for name in keys(data_dfs)]...;
    overwrite,
  )
end

function loadModel(filename::String)
  read_from_file(filename)
end

