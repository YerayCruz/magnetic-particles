using DataStructures

info = DefaultDict(0.)
function file_reader(file)

    f = open(file, "r")
    regex = r"(\w+) = (\d+)"

    s = read(f, String)
    matches = eachmatch(regex, s)
    for match in matches
        info["$(match.captures[1])"] = parse(Float64, match.captures[2])
        println(info)
    end
end

file_reader("./src/file_reader/textfile.txt")
