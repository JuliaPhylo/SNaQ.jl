# Search log files.

"""
Helper function to log the message `msg` to file `logfile`.
"""
function logtext(logfile::String, msg::String)
    logfile == "" && return
    # get the current time and format it
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    # open the file in append mode and write the log line
    open(logfile, "a") do io
        println(io, "[$timestamp] $msg")
    end
end


"""
Helper function to log proposed and accepted moves to `logfile`.
"""
function logmoves(logfile::String, moves_prop::Dict, moves_acc::Dict, moves_PL::Dict)
    all_keys = sort(collect(keys(moves_prop)))
    min_width::Int = maximum([max(length(string(k)), length(string(moves_prop[k])), length(string(moves_acc[k]))) for k in all_keys])+2
    function expand(str::String)::String
        ret::String = str
        for j = 1:(min_width - length(str))
            ret *= " "
        end
        return ret
    end

    msg::String = repeat("-", 12) * "MOVE ACCEPTANCE RATES" * repeat("-", 12) * "\n"
    msg *= expand("\tmove:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(k))
    end
    msg *= " |\n" * expand("\tproposals:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(moves_prop[k]))
    end
    msg *= " |\n" * expand("\taccepted:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(moves_acc[k]))
    end
    msg *= " |\n" * expand("\taccept %:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(round(100 * moves_acc[k] / moves_prop[k], digits=0)))
    end
    msg *= " |\n" * expand("\tmean SUCC ΔPL:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(round(mean(moves_PL[k]), digits=2)))
    end
    msg *= " |\n" * expand("\tmax SUCC ΔPL:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(round(maximum(moves_PL[k], init=0), digits=2)))
    end
    msg *= " |\n" * expand("\tmean ALL ΔPL:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(round(sum(moves_PL[k]) / moves_prop[k], digits=2)))
    end
    msg *= " |\n\n"

    logtext(logfile, msg)
end


logmessage(filename::String, msg::String) = remotecall_fetch(writelogmessage, 1, filename, msg)


"""
    @logmessage filename msg

[`logmessage`](@ref), but the message expression is only evaluated when there is a file
to write it to. Several log messages interpolate `writenewick(...)` of the whole network,
which is expensive on a large one -- and was paid on every `search` call even with
logging turned off, since the string is built before `logmessage` can discard it.
"""
macro logmessage(filename, msg)
    return :(local f = $(esc(filename)); f == "" || logmessage(f, $(esc(msg))))
end


currenttime() = Dates.format(now(), "HH:MM:SS yyyy-mm-dd")


function writelogmessage(filename::String, msg::String)
    filename == "" && return
    open("$(filename).log", "a+") do f
        println(f, msg)
    end
end


function timeelapsed(totaltime::Float64)::String
    seconds::Int = Int(round(totaltime))
    minutes::Int = (seconds ÷ 60) % 60
    hours::Int   = (seconds ÷ 3600) % 24
    days::Int    = seconds ÷ 86400
    seconds      = seconds % 60
    if days > 0
        return "$days days, $hours hours, $minutes minutes and $seconds seconds"
    elseif hours > 0
        return "$hours hours, $minutes minutes and $seconds seconds"
    elseif minutes > 0
        return "$minutes minutes and $seconds seconds"
    else
        return "$seconds seconds"
    end
end
