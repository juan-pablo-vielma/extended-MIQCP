
function userTime()
	return  ccall( (:user_time, "rusagetimer"), Cdouble, ())
end

function userAndSystemTime()
	return  ccall( (:user_and_system_time, "rusagetimer"), Cdouble, ())
end

