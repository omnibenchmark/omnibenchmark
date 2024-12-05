-- Python_3.12.6_Clustering

help([[
    Python 3.12.6 environment with clustering libraries
    Includes scikit-learn, pandas, scipy, and numpy
]])

whatis("Name: Python_3.12.6_Clustering")
whatis("Version: 1.0.0")
whatis("Category: OmniBenchmark")
whatis("Description: Python environment for omnibenchmark clustering")

-- Python Configuration
setenv("PYTHON_VERSION", "3.12.6")

-- Define the Python environment path (adjust to where Python 3.12.6 is installed)
local python_version = "3.12"  -- Modify this to the actual Python path
local python_bin = capture("python" .. python_version .. " -c 'import sys; print(sys.executable)'")
local python_home = capture("python" .. python_version .. " -c 'import sys; print(sys.prefix)'")
local python_lib = python_home .. "/lib/python" .. python_version
local python_site_packages = python_lib .. "/site-packages"

-- Set up the environment to use Python 3.12.6
prepend_path("PATH", python_bin)
prepend_path("LD_LIBRARY_PATH", python_lib)
prepend_path("PYTHONPATH", python_site_packages)

-- Load necessary Python libraries (assuming they are installed in the specified path)
local required_libraries = {
    "scikit-learn",
    "pandas",
    "scipy",
    "numpy"
}

-- Make sure the required Python packages are accessible by adding their site-packages path
for _, library in ipairs(required_libraries) do
    local library_path = "/usr/lib/python" .. python_version .. "/site-packages/" .. library
    prepend_path("PYTHONPATH", library_path)
end


-- Print information about the environment
LmodMessage("Python 3.12.6 and clustering libraries have been successfully loaded.")

---
-- Function to retrieve console output
--
function capture(cmd, raw)
    local handle = assert(io.popen(cmd, 'r'))
    local output = assert(handle:read('*a'))

    handle:close()

    if raw then
        return output
    end

    output = string.gsub(
        string.gsub(
            string.gsub(output, '^%s+', ''),
            '%s+$',
            ''
        ),
        '[\n\r]+',
        ' '
    )

   return output
end