$ErrorActionPreference = "Stop"
$root = "C:\Users\53461\OneDrive\Documentos\Patch-seq_Mouse_Acb_MSN_T-type_Visualization"
Set-Location -LiteralPath $root
python "scripts\build_macaque_M4_frozen_input_audit.py"
python "scripts\validate_macaque_M4_full_ml500.py" --repeats 500 --jobs -1
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "scripts\reproduce_macaque_M4_frozen_gc.R"
