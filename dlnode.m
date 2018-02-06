% Object-Oriented Programming 
% Doubly-Linked List 
%
% Example:
% Makes doubly-linked list with three nodes with data values 1, 2, and 3:
% n1 = dlnode(1);
% n2 = dlnode(2);
% n3 = dlnode(3); 
% n2.insertAfter(n1)    insert n2 after n1
% n3.insertAfter(n2)    insert n3 after n2
% n1.Next               points to n2
% n2.Next.Prev          points back to n2
% n1.Next.Next          points to n3
% n3.Prev.Prev          points to n1 

classdef dlnode < handle
% DLNODE A class to represent a doubly-linked list. 
% Multiple dlnode objects may be linked together to create linked lists
% Each node contains a piece of data and provides access to the next and
% previous nodes. 
    properties 
        Data
    end
    properties(SetAccess = private)
        Next
        Prev
    end
    
    methods
        function node = dlnode(Data)
        % DLNODE Constructs a dlnode object.
        if nargin > 0 
            node.Data = Data; 
        end
        end
        
        function insertAfter(newNode, nodeBefore)
            % insertAfter Inserts newNode after newBefore
            disconnect(newNode);
            newNode.Next = nodeBefore.Next;
            newNode.Prev = nodeBefore;
            if ~isempty(nodeBefore.Next)
            nodeBefore.Next.Prev = newNode;
            end
            nodeBefore.Next = newNode; 
        end
        
        function insertBefore(newNode, nodeAfter)
            % insertBefore Inserts newNode before nodeAfter 
            disconnect(newNode);
            newNode.Next = nodeAfter;
            newNode.Prev = nodeAfter.Prev;
            if ~isempty(nodeAfter.Prev)
                nodeAfter.Prev.Next = newNode;
            end
            nodeAfter.Prev = newNode;
        end
        
        function disconnect(node)
            % DISCONNECT Removes a node from a linked list
            % The node can be reconnected or moved to a different list
            if ~isempty(node.Prev)
                node.Prev.Next = node.Next;
            end
            if ~isempty(node.Next)
                node.Next.Prev = node.Prev;
            end
            node.Next = [];
            node.Prev = [];
        end
        
        function delete(node)
            % DELETE Deletes a dlnode from a linked list
            disconnect(node);
        end
        
        function disp(node)
            % DISP Displays a link node
            disp('Doubly-linked list node with data:');
            disp(node.Data);
        end
    end
end


    
    
    
            
    