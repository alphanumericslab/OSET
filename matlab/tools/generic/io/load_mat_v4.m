function [data, name, T] = load_mat_v4(filename)
    % LOAD_MATV4 Read and extract data from a MAT-4 file.
    %
    %   [DATA, NAME] = load_mat_v4(FILENAME) reads a MAT-4 file specified by
    %   FILENAME and extracts the data stored in the file. It returns the
    %   extracted matrix data in DATA and the name of the matrix in NAME.
    %
    %   Example:
    %   [data, name] = load_mat_v4('0284_001_004_ECG.mat');
    %
    %   Reference:
    %   The function is based on the Level 4 MAT-File Matrix Header Format.
    %   For more details, refer to pages 1-25 throgh 1-28:
    %   https://www.mathworks.com/help/pdf_doc/matlab/matfile_format.pdf
    %
    %   Reza Sameni, June 2023
    %

    fid = fopen(filename, 'rb');

    if fid == -1
        error('File could not be opened.');
    end

    try
        % Read the matrix header
        header = fread(fid, 5, 'uint32');

        % Extract the relevant header fields
        type = header(1);
        M = floor(mod(type/1000,10));  % M indicates the numeric format of binary numbers
        % O = floor(mod(type/100,10));   % O is always 0 (zero) and is reserved for future use
        P = floor(mod(type/10,10));    % P indicates the format of the data stored
        T = floor(mod(type,10));       % T indicates the matrix type

        mrows = header(2);     % The row dimension contains the number of rows in the matrix
        ncols = header(3);     % The column dimension contains the number of columns in the matrix
        imagf = header(4);     % The imaginary flag: 1 if matrix has imaginary part, 0 if only real data
        namelen = header(5);   % The name length contains an integer with 1 plus the length of the matrix name

        % Read the matrix name
        name = fread(fid, namelen - 1, 'uint8=>char')';

        % Determine the byte order based on M field
        if M == 0
            byteOrder = 'l';   % IEEE Little Endian
        elseif M == 1
            byteOrder = 'b';   % IEEE Big Endian
        else
            error('Unsupported byte order.');
        end

        % Read the real part of the matrix
        real_data = fread(fid, mrows * ncols, getDataType(P), 0, byteOrder);

        % Read the imaginary part if present
        if imagf ~= 0
            imag_data = fread(fid, mrows * ncols, getDataType(P), 0, byteOrder);
            data = complex(reshape(real_data, mrows, ncols), reshape(imag_data, mrows, ncols));
        else
            data = reshape(real_data, mrows, ncols);
        end

        % Close the file
        fclose(fid);
    catch ex
        % Close the file in case of an error
        fclose(fid);
        rethrow(ex);
    end
end

function dtype = getDataType(type)
    % getDataType maps the P value to the corresponding data type
    switch type
        case 0
            dtype = 'float64';   % double-precision (64-bit) floating-point numbers
        case 1
            dtype = 'float32';   % single-precision (32-bit) floating-point numbers
        case 2
            dtype = 'int32';     % 32-bit signed integers
        case 3
            dtype = 'int16';     % 16-bit signed integers
        case 4
            dtype = 'uint16';    % 16-bit unsigned integers
        case 5
            dtype = 'uint8';     % 8-bit unsigned integers
        otherwise
            error('Unsupported data type.');
    end
end
