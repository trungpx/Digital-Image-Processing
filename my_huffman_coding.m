
%% How to run: Just click "Run" button and hit "Select folder"

function my_huffman_coding()
clc;
clear;
close all;
imgdir = uigetdir('Test_images');
%% Loading directory, gray image
file = fopen(fullfile(imgdir,'\bridge_gray_256x256.raw'),'rb');
gray_image = fread(file,fliplr([256,256]),'*uint8')';
fclose(file); gray_image = double(gray_image);

%% starting
disp('Huffman coding');
disp('Computing the hufman code...');

% Computes Hufmancode
[Huffmancode,p,H,CL,eff] = my_huffman_func(gray_image);
Huffmanword = Huffman_words(p);

% Displays figures
figure; imshow(gray_image,[]); title('Problem 2) Origional'); % show original image
disp('Huffman words table: ');
T = cell2table(Huffmanword);
Tnew =array2table(0:255);
Tnew.Properties.VariableNames = T.Properties.VariableNames;
[Tnew; T]
%disp(['Probabilities: ',num2str(p)]);
%disp(['Huffman_code_sequence: ',num2str(Huffmancode)]);
disp(['Entropy of original image: ',num2str(H)]);
disp(['Average code length: ',num2str(CL)]);
disp(['Compress efficiency: ',num2str(eff)]);

end

%% huffmancode function
function [hcode, probability, H, codelength, efficiency] = my_huffman_func(I)
%H is entropy
%size of the image
[m,n]=size(I);
Totalcount=m*n;
%variables using to find the probability
cnt=1;
sigma=0;
%Symbols for an image
symbols = 0:255;
%computing the cumulative probability.
for i=0:255
    k=I==i;
    count(cnt)=sum(k(:));
    %pro array is having the probabilities
    pro(cnt)=count(cnt)/Totalcount;
    sigma=sigma+pro(cnt);
    cnt=cnt+1;
end
% keep aligned the vector of symbols and the vector of their probabilities
ind=find(pro>0);
ExistingSymbols = symbols(ind);
probability = pro(ind);
%Huffman code Dictionary
dict = my_huffmandict(ExistingSymbols,probability);
% Entropy
H = -sum(probability.*log2(probability));
% Check the theoretical average codelength
CL_symb = zeros(size(ExistingSymbols));
for i =1:length(ExistingSymbols)
    CL_symb(i) = length(dict{i,2});
end
codelength = sum(probability.*CL_symb);
% find hcode
MCL = max(CL_symb);
N=length(ExistingSymbols);
CodeExt = -1*ones(MCL,N);
EncodedMatrix = zeros(MCL,length(I(:)));
for i =1:N
    CL_s  = length(dict{i,2});
    CodeExt(1:CL_s,i) = dict{i,2}';
end
for i =1:N
    ind = find(I(:) == ExistingSymbols(i));
    EncodedMatrix(1:MCL,ind) = CodeExt(1:MCL,i)*ones(1,length(ind));
end
ind = find( EncodedMatrix(:) ~= -1);

hcode = EncodedMatrix(ind)';
efficiency = H/codelength;
end

%% function for new my_huffmandict
function [dict] = my_huffmandict(sig, prob)
n_ary = [];
variance = 'max';
if isempty(n_ary)
    n_ary = 2; % default value is binary encryption
end
% Make sure that internally all vectors are represented as column vectors
m = size(sig);
if( m(1) == 1 )
    sig = sig';
end
prob = prob(:);
% Make sure that the input symbols are in a cell array format
if ~iscell(sig)
    sig = num2cell(sig) ;
end
% === Create the dictionary ===
% Create tree nodes with the signals and the corresponding probabilities
huff_tree = struct('signal', [], 'probability', [],...
    'child', [], 'code', [], 'origOrder', -1);
for i=1:length( sig )
    huff_tree(i).signal = sig{i};
    huff_tree(i).probability = prob(i);
    huff_tree(i).origOrder = i;
end
% Sort the signal and probability vectors based on ascending order of
% probability
[~, i] = sort(prob);
huff_tree = huff_tree(i);
huff_tree = create_tree(huff_tree, n_ary, variance); % create a Huffman tree
[~,dict,~] = create_dict(huff_tree,{},0); % create the codebook
% The next few lines of code are to sort the dictionary.
% If sorting based on original order then use dict{:,4}.
[~,dictsortorder] = sort([dict{:,4}]);
lenDict = length(dictsortorder);
finaldict = cell(lenDict, 2);
for i=1:length(dictsortorder)
    finaldict{i,1} = dict{dictsortorder(i), 1};
    finaldict{i,2} = dict{dictsortorder(i), 2};
end
dict = finaldict;
end
%% create tree
function huff_tree = create_tree(huff_tree, n_ary, variance)
% if the length of huff_tree is 1, it implies there is no more than one
% node in the array of nodes. This is the termination condition for the
% recursive loop
numRemNodes = length(huff_tree);
if( numRemNodes <= 1)
    return;
end
temp = struct('signal', [], 'probability', 0, ...
    'child', [], 'code', []);
numNodesToComb = rem(numRemNodes-1, n_ary-1) + 1;
if numNodesToComb == 1 % Must be true except for the first round
    numNodesToComb = n_ary;
end

for i = 1:numNodesToComb
    if isempty(huff_tree), break; end
    temp.probability = temp.probability + huff_tree(1).probability; % for ascending order
    temp.child{i} = huff_tree(1);
    temp.origOrder = -1;
    huff_tree(1) = [];
end
if( strcmpi(variance, 'min') == 1 )
    huff_tree = insertMinVar(huff_tree, temp);
else
    huff_tree = insertMaxVar(huff_tree, temp);
end
% create a Huffman tree from the reduced number of free nodes
huff_tree = create_tree(huff_tree, n_ary, variance);
return

end
%% insert max value
function huff_tree = insertMaxVar(huff_tree, newNode)
i = 1;
while i <= length(huff_tree) && ...
        newNode.probability > huff_tree(i).probability
    i = i+1;
end
huff_tree = [huff_tree(1:i-1) newNode huff_tree(i:end)];
end
%% insert min value
function huff_tree = insertMinVar(huff_tree, newNode)
i = 1;
while i <= length(huff_tree) && ...
        newNode.probability >= huff_tree(i).probability
    i = i+1;
end
huff_tree = [huff_tree(1:i-1) newNode huff_tree(i:end)];
end
%% create dict
function [huff_tree,dict,total_wted_len] = create_dict(huff_tree,dict,total_wted_len)
% Check if the current node is a leaf node If it is, then add the signal on
% this node and its corresponding code to the dictionary global n_ary
if isempty(huff_tree.child)
    dict{end+1,1} = huff_tree.signal;
    dict{end, 2} = huff_tree.code;
    dict{end, 3} = length(huff_tree.code);
    dict{end, 4} = huff_tree.origOrder;
    total_wted_len = total_wted_len + length(huff_tree.code)*huff_tree.probability;
    return;
end
num_childrens = length(huff_tree.child);
for i = 1:num_childrens
    huff_tree.child{i}.code = [huff_tree(end).code, (num_childrens-i)];
    [huff_tree.child{i}, dict, total_wted_len] = ...
        create_dict(huff_tree.child{i}, dict, total_wted_len);
end
end
%% corewords
function [h]=Huffman_words(p)
% Huffman code generator gives a Huffman code matrix h, 
%  average codeword length L & entropy H
% for a source with probability vector p given as argin(1) 
zero_one=['0'; '1']; 
M=length(p);  N=M-1; p=p(:); % Make p a column vector
h={zero_one(1),zero_one(2)};
if M>2
  pp(:,1)=p;
  for n=1:N
     % To sort in descending order
     [pp(1:M-n+1,n),o(1:M-n+1,n)]=sort(pp(1:M-n+1,n),1,'descend'); 
     if n==1, ord0=o; end  % Original descending order
     if M-n>1, pp(1:M-n,n+1)=[pp(1:M-1-n,n); sum(pp(M-n:M-n+1,n))]; end
  end
  for n=N:-1:2
     tmp=N-n+2; oi=o(1:tmp,n);
     for i=1:tmp, h1{oi(i)}=h{i}; end
     h=h1;   h{tmp+1}=h{tmp};
     h{tmp}=[h{tmp} zero_one(1)]; 
     h{tmp+1}=[h{tmp+1} zero_one(2)];
  end
  for i=1:length(ord0), h1{ord0(i)}=h{i}; end
  h=h1;
end
end