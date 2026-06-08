# This function helps safely redirect stdout and stderr temporarily
function safely_redirect_output(f::Function)
    redirect_stdout(devnull) do 
        redirect_stderr(devnull) do 
            f()
        end
    end
end