# export/xml.jl -- a deliberately tiny XML well-formedness checker, used to validate our
# own SVG and 3MF output without pulling in a parser. It checks: a single root
# element, balanced and correctly nested tags, quoted attribute values, no stray
# '<' or '>' in text, and well-formed comments / declarations / CDATA.
#
# Works on code units (bytes) with C-locale (ASCII) character classes; byte positions
# in messages are 0-based.

_xml_isalpha(ch::UInt8) = (UInt8('a') <= ch <= UInt8('z')) || (UInt8('A') <= ch <= UInt8('Z'))
_xml_isdigit(ch::UInt8) = UInt8('0') <= ch <= UInt8('9')
_xml_isspace(ch::UInt8) = ch == 0x20 || (0x09 <= ch <= 0x0d)
_name_start(ch::UInt8) = _xml_isalpha(ch) || ch == UInt8('_') || ch == UInt8(':')
_name_char(ch::UInt8) = _xml_isalpha(ch) || _xml_isdigit(ch) || ch == UInt8('_') ||
                        ch == UInt8(':') || ch == UInt8('-') || ch == UInt8('.')

# position of `pat` in `s` at or after byte index i (1-based), or 0
function _bfind(s::AbstractString, pat::AbstractString, i::Int)
    i > ncodeunits(s) && return 0
    r = findnext(pat, s, i)
    return r === nothing ? 0 : first(r)
end
_bfind(s::AbstractString, ch::Char, i::Int) = _bfind(s, string(ch), i)
_startswith_at(s, i, pat) = i + ncodeunits(pat) - 1 <= ncodeunits(s) &&
                            SubString(s, i, i + ncodeunits(pat) - 1) == pat

"""
    xml_well_formed(text) -> (ok, err)

True iff `text` is a well-formed XML document by the rules above.
On failure `err` receives a human-readable reason.
"""
function xml_well_formed(s::AbstractString)
    n = ncodeunits(s)
    stack = String[]
    roots = 0
    i = 1
    fail(msg, pos) = (false, msg * " at byte " * string(pos - 1))
    at(k) = codeunit(s, k)
    while i <= n
        if at(i) != UInt8('<')
            at(i) == UInt8('>') && return fail("bare '>' in character data", i)
            if at(i) == UInt8('&')
                semi = _bfind(s, ';', i)
                (semi == 0 || semi - i > 12) && return fail("unterminated entity reference", i)
                i = semi + 1
                continue
            end
            (isempty(stack) && !_xml_isspace(at(i))) &&
                return fail("character data outside the root element", i)
            i += 1
            continue
        end
        # markup
        if _startswith_at(s, i, "<!--")
            e = _bfind(s, "-->", i + 4)
            e == 0 && return fail("unterminated comment", i)
            i = e + 3
            continue
        end
        if _startswith_at(s, i, "<![CDATA[")
            e = _bfind(s, "]]>", i + 9)
            e == 0 && return fail("unterminated CDATA section", i)
            i = e + 3
            continue
        end
        if _startswith_at(s, i, "<?")
            e = _bfind(s, "?>", i + 2)
            e == 0 && return fail("unterminated processing instruction", i)
            isempty(stack) || return fail("processing instruction inside an element", i)
            i = e + 2
            continue
        end
        if _startswith_at(s, i, "<!")  # DOCTYPE and friends
            e = _bfind(s, '>', i + 2)
            e == 0 && return fail("unterminated declaration", i)
            i = e + 1
            continue
        end
        closing = (i + 1 <= n && at(i + 1) == UInt8('/'))
        j = i + (closing ? 2 : 1)
        (j > n || !_name_start(at(j))) && return fail("malformed tag name", i)
        name_beg = j
        while j <= n && _name_char(at(j))
            j += 1
        end
        tag = String(SubString(s, name_beg, j - 1))

        if closing
            while j <= n && _xml_isspace(at(j))
                j += 1
            end
            (j > n || at(j) != UInt8('>')) && return fail("malformed end tag", i)
            isempty(stack) && return fail("end tag </$tag> with no open element", i)
            stack[end] != tag && return fail("end tag </$tag> closes <$(stack[end])>", i)
            pop!(stack)
            i = j + 1
            continue
        end
        # attributes
        self_closing = false
        while true
            while j <= n && _xml_isspace(at(j))
                j += 1
            end
            j > n && return fail("unterminated start tag", i)
            if at(j) == UInt8('>')
                j += 1
                break
            end
            if at(j) == UInt8('/')
                (j + 1 > n || at(j + 1) != UInt8('>')) && return fail("malformed self-closing tag", j)
                self_closing = true
                j += 2
                break
            end
            _name_start(at(j)) || return fail("malformed attribute name", j)
            while j <= n && _name_char(at(j))
                j += 1
            end
            while j <= n && _xml_isspace(at(j))
                j += 1
            end
            (j > n || at(j) != UInt8('=')) && return fail("attribute without a value", j)
            j += 1
            while j <= n && _xml_isspace(at(j))
                j += 1
            end
            (j > n || (at(j) != UInt8('"') && at(j) != UInt8('\''))) &&
                return fail("unquoted attribute value", j)
            quote_ch = Char(at(j))
            j += 1
            e = _bfind(s, quote_ch, j)
            e == 0 && return fail("unterminated attribute value", j)
            lt = _bfind(s, '<', j)
            (lt != 0 && lt < e) && return fail("'<' inside an attribute value", j)
            j = e + 1
        end
        if !self_closing
            isempty(stack) && (roots += 1)
            push!(stack, tag)
        elseif isempty(stack)
            roots += 1
        end
        i = j
    end
    isempty(stack) || return fail("unclosed element <$(stack[end])>", n + 1)
    roots != 1 && return (false, "document has $roots root elements, expected 1")
    return true, ""
end

"""
Number of occurrences of the element `name` (counting `<name ...>` and
`<name/>`, not closing tags). Returns 0 when absent.
"""
function xml_count_elements(s::AbstractString, name::AbstractString)
    count_ = 0
    open_ = "<" * name
    n = ncodeunits(s)
    i = 1
    while (i = _bfind(s, open_, i)) != 0
        j = i + ncodeunits(open_)
        if j <= n && (_xml_isspace(codeunit(s, j)) || codeunit(s, j) == UInt8('>') ||
                      codeunit(s, j) == UInt8('/'))
            count_ += 1
        end
        i = j
    end
    return count_
end

"""The values of attribute `attr` on every `<name ...>` start tag, in order."""
function xml_attribute_values(s::AbstractString, name::AbstractString, attr::AbstractString)
    out = String[]
    open_ = "<" * name
    n = ncodeunits(s)
    i = 1
    while (i = _bfind(s, open_, i)) != 0
        j = i + ncodeunits(open_)
        if j > n || !(_xml_isspace(codeunit(s, j)) || codeunit(s, j) == UInt8('>') ||
                      codeunit(s, j) == UInt8('/'))
            i = j
            continue
        end
        e = _bfind(s, '>', j)
        e == 0 && break
        tag = SubString(s, j, e - 1)
        a = _bfind(tag, attr * "=\"", 1)
        if a != 0
            vb = a + ncodeunits(attr) + 2
            ve = _bfind(tag, '"', vb)
            ve != 0 && push!(out, String(SubString(tag, vb, ve - 1)))
        end
        i = e + 1
    end
    return out
end
