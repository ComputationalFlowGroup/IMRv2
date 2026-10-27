% file f_format_lines.m
% brief contains function f_format_lines

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function re-indents one file's lines, given as a cell array
% of char rows (no trailing newline), and returns the re-indented lines
% the same way. It is the core of s_format.m, factored out so that
% s_lint.m can run it on a fixture and catch a regression that mangles
% code rather than only re-indenting it: it must never change the number
% of lines, and it must never alter what a line says, only its leading
% whitespace. In particular a quoted string or a comment is opaque to it:
% it is never split, searched for keywords, or counted for semicolons,
% because telling a statement-separating ';' from one inside a string or
% a comment needs a real tokenizer, not a character count.
function formattedLines = f_format_lines(lines)

    indentLevel = 0;
    formattedLines = cell(numel(lines),1);
    count = 1;
    indentStep = '    ';
    multiLineContinuation = false;

    for j = 1:numel(lines)
        originalLine = lines{j};
        line = strtrim(originalLine);
        % Preserve blank lines instead of adding trailing indentation.
        if isempty(line)
            formattedLines{count} = originalLine;
            count = count + 1;
            continue;
        end

        if startsWith(line, "end")
            indentLevel = max(indentLevel - 1, 0);
        end

        if multiLineContinuation
            formattedLines{count} = [repmat(indentStep, 1, indentLevel + 1), line];
        else
            if any(startsWith(line, ["else", "catch"]))
                formattedLines{count} = [repmat(indentStep, 1, max(indentLevel - 1, 0)), line];
            else
                formattedLines{count} = [repmat(indentStep, 1, indentLevel), line];
            end
        end
        count = count + 1;

        multiLineContinuation = endsWith(line, "...");

        if ~multiLineContinuation && any(startsWith(line, ["if ", "parfor ", "for ", "while ", "switch ", "function", "try"]))
            indentLevel = indentLevel + 1;
        end
    end

    formattedLines = formattedLines(1:count-1);
    lastNonEmptyIdx = find(~cellfun(@isempty, formattedLines), 1, 'last');
    formattedLines = formattedLines(1:lastNonEmptyIdx);

end
