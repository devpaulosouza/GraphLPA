#include <iostream>
#include <stdio.h>

#define WHITE (0)
#define GREY  (1)
#define BLACK (2)
#define INT32_MAX		(2147483647)
#define nullptr NULL
#define IOFB "IndexOutOfBoundsException"

using namespace std;


class Vertex {
    int label;
    int nEdges;
    Vertex *edges;
    int inDegree, outDegree;
    // BFS vars
    short color;
    int bfsDistance;
    Vertex *dad;
public:
    Vertex();

    explicit Vertex(int label);

    int getLabel();

    Vertex getVertex(int label);

    int getNEdges();

    Vertex *edgeAt(int position);

    void addEdge(Vertex next);

    void removeEdge(int label);

    void printEdges();

    short getColor() const;

    void setColor(short color);

    int getBfsDistance() const;

    void setBfsDistance(int bfsDistance);

    Vertex *getDad() const;

    void setDad(Vertex *dad);

    void setLabel(int label);

    int getInDegree() const;

    void incrementInDregree();

    int getOutDegree() const;

    void incrementOutDegree();
};

Vertex::Vertex() {
    this->label = -1;
    this->nEdges = 0;
    this->edges = new Vertex[0]; // NOLINT
    this->color = WHITE;
    this->bfsDistance = INT32_MAX;
    this->dad = nullptr;
    this->inDegree = this->outDegree = 0;
}

Vertex::Vertex(int label) {
    this->edges = new Vertex[0]; // NOLINT
    this->label = label;
    this->nEdges = 0;
    this->bfsDistance = INT32_MAX;
    this->color = WHITE;
    this->dad = nullptr;
    this->inDegree = this->outDegree = 0;
}

int Vertex::getLabel() {
    return this->label;
}

Vertex Vertex::getVertex(int label) {
    return edges[label];
}

int Vertex::getInDegree() const {
    return inDegree;
}

int Vertex::getOutDegree() const {
    return outDegree;
}

void Vertex::addEdge(Vertex next) {
    bool vertexAlreadyExists = false;
    for (int i = 0; i < this->nEdges; ++i) {
        vertexAlreadyExists |= this->edges[i].getLabel() == next.getLabel();
    }
    if (vertexAlreadyExists) {
        //cout << "Vertex " << this->label << " - " << next.label << " already exists" << endl;
        return;
    }

    Vertex *edges = new Vertex[++this->nEdges]; // NOLINT
    for (int i = 0; i < this->nEdges - 1; ++i) {
        edges[i] = this->edges[i];
    }
    this->edges = new Vertex[nEdges];
    for (int i = 0; i < this->nEdges - 1; ++i) {
        this->edges[i] = edges[i];
    }
    this->edges[this->nEdges - 1] = next;
}

void Vertex::printEdges() {
    if (nEdges == 0) {
        cout << "No edges." << endl;
    } else {
        cout << "Edges: ";
        for (int i = 0; i < nEdges; ++i) {
            // pretty-fy ^^
            if (i != nEdges - 1) {
                cout << "(" << this->label << ", " << this->edges[i].getLabel() << "), ";
            } else {
                cout << "(" << this->label << ", " << this->edges[i].getLabel() << ")";
            }
        }
        cout << endl;
    }
}

void Vertex::removeEdge(int label) {
    int i;
    for (i = 0; i < this->nEdges && this->edges[i].label != label; ++i) {};
    if (i < nEdges) {
        Vertex *newEdges = new Vertex[--this->nEdges]; // NOLINT
        for (int j = 0; j < this->nEdges; ++j) {
            if (j != i) {
                newEdges[j] = this->edges[j];
            }
        }
    } else {
        cout << "Edge not found. From: " << this->label << " to " << label << endl;
    }
}

short Vertex::getColor() const {
    return color;
}

void Vertex::setColor(short color) {
    Vertex::color = color;
}

int Vertex::getBfsDistance() const {
    return bfsDistance;
}

void Vertex::setBfsDistance(int bfsDistance) {
    Vertex::bfsDistance = bfsDistance;
}

Vertex *Vertex::getDad() const {
    return dad;
}

void Vertex::setDad(Vertex *dad) {
    Vertex::dad = dad;
}

int Vertex::getNEdges() {
    return nEdges;
}

Vertex *Vertex::edgeAt(int position) {
    if (position >= this->nEdges) {
        cout << "Edge not found at position " << position << endl;
        cout << "Range expected 0.." << this->nEdges << endl;
        __throw_logic_error("EdgeNotFoundException");
    }
    return &this->edges[position];
}

void Vertex::setLabel(int label) {
    Vertex::label = label;
}

void Vertex::incrementInDregree() {
    this->inDegree++;
}

void Vertex::incrementOutDegree() {
    this->outDegree++;
}

class Queue {
    Vertex *queue;
    int size;
public:
    Queue();

    Vertex *unqueue();

    void enqueue(Vertex vertex);

    int getSize() const;
};

void Queue::enqueue(Vertex vertex) {
    Vertex *aux = new Vertex[++this->size];

    for (int i = 0; i < this->size - 1; ++i) {
        aux[i] = queue[i];
    }
    aux[this->size - 1] = vertex;

    queue = aux;
    /*
    Vertex *aux = new Vertex[vertex->getNEdges()];
    int i;
    for (i = 0; i < quantity; ++i) {
    if (vertex->edgeAt(i).getColor() == WHITE) {
    aux[i] = vertex->edgeAt(i);
    }
    }

    Vertex *newQueue = new Vertex[quantity + i];
    for (int j = 0; j < this->size; ++j) {
    newQueue[j] = this->queue[j];
    }
    for (int j = size; j < size+quantity; ++j) {
    newQueue[j] = aux[j-size];
    }

    this->size += i;
    queue = newQueue;
    */
}

int Queue::getSize() const {
    return this->size;
}

Queue::Queue() {
    this->queue = new Vertex[0];
    this->size = 0;
}

Vertex *Queue::unqueue() {
    Vertex *newQueue = new Vertex[--this->size];
    Vertex *callback;
    if (this->size == -1) {
        __throw_length_error("Queue is empty");
    }
    callback = &queue[0];
    for (int i = 0; i < this->size; ++i) {
        newQueue[i] = queue[i + 1];
    }
    queue = newQueue;
    return callback;
}

class Graph {
    int size;
    Vertex *allVertex;
// BFS vars
    Vertex *queue;
public:

    explicit Graph(int nVertex);

    void print();

    void addEdge(int predecessor, int successor);

    void removeEdge(int predecessor, int successor);

    bool bfs(int begin);

    bool irEVir ();
};

Graph::Graph(int nVertex) {
    Vertex *allVertex = new Vertex[nVertex]; // NOLINT
    for (int i = 0; i < nVertex; ++i) {
        allVertex[i] = Vertex(i);
    }
    this->allVertex = allVertex;
    size = nVertex;
    queue = nullptr;
}

void Graph::print() {
    for (int i = 0; i < this->size; ++i) {
        cout << "Vertex: " << this->allVertex[i].getLabel() << endl;
        cout << "InDegree: " << this->allVertex[i].getInDegree() << endl;
        cout << "OutDegree: " << this->allVertex[i].getOutDegree() << endl;
        cout << "Color:  " << this->allVertex[i].getColor() << endl;
        cout << "bfs:    " << this->allVertex[i].getBfsDistance() << endl;
        this->allVertex[i].printEdges();
        cout << endl;
    }
}

void Graph::addEdge(int predecessor, int successor) {
    if (predecessor >= size || successor >= size) {
        cout << "Index out of bounds in method addEdge." << endl;
        cout << "predecessor = " << predecessor << ", successor = " << successor << endl;
        cout << "Range supported: 0.." << size - 1 << endl;
        __throw_out_of_range(IOFB);
    }
    allVertex[predecessor].addEdge(allVertex[successor]);
    allVertex[predecessor].incrementOutDegree();
    allVertex[successor].incrementInDregree();
}

void Graph::removeEdge(int predecessor, int successor) {
    if (predecessor >= size || successor >= size) {
        cout << "Index out of bounds in method removeEdge." << endl;
        cout << "predecessor = " << predecessor << ", successor = " << successor << endl;
        cout << "Range supported: 0.." << size - 1 << endl;
        __throw_out_of_range(IOFB);
    }
    allVertex[predecessor].removeEdge(successor);
}


bool Graph::bfs(int begin) {
    bool res = true;
    if (begin >= this->size) {
        cout << "Index out of bounds in method bfs." << endl;
        cout << "begin = " << begin << endl;
        cout << "Range supported: 0.." << size - 1 << endl;
        __throw_out_of_range(IOFB);
    }
    Vertex *root = &allVertex[begin];

    // reset values
    for (int i = 0; i < this->size; ++i) {
        this->allVertex[i].setBfsDistance(INT32_MAX);
        this->allVertex[i].setColor(WHITE);
        this->allVertex[i].setDad(nullptr);
    }

    // mark as visited
    root->setColor(GREY);
    root->setBfsDistance(0);
    root->setDad(nullptr);
    // create empty queue
    Queue queue;
    // enqueue root
    queue.enqueue(*root);

    while (queue.getSize() != 0) {
        // get first vertex
        // Semente dos deuses para o algoritmo funcionar
        Vertex *actual = &allVertex[queue.unqueue()->getLabel()];
        Vertex *neightbor;
        for (int i = 0; i < actual->getNEdges(); ++i) {
            // Semente dos deuses para o algoritmo funcionar
            neightbor = &allVertex[actual->edgeAt(i)->getLabel()];
            if (neightbor->getColor() == WHITE) {
                neightbor->setColor(GREY);
                neightbor->setBfsDistance(actual->getBfsDistance() + 1);
                neightbor->setDad(actual);
                queue.enqueue(*neightbor);
            }
        }
        actual->setColor(BLACK);
    }

    for (int i = 0; i < this->size && res; ++i) {
        res &= allVertex[i].getColor()!=WHITE;
    }
    return res;
}

bool Graph::irEVir() {
    bool res = true;
    for (int i = 0; i < this->size && res; ++i) {
        res &= bfs(i);
    }
    return res;
}



int main() {
    int streets, corners;
    int v1, v2, ways;
    Graph *city;

    cin >> corners >> streets;

    while (corners != 0 && streets != 0) {
        city = new Graph(corners);

        for (int i = 0; i < streets; ++i) {
            cin >> v1 >> v2 >> ways;

            v1--; v2--;
            if (ways == 2) {
                city->addEdge(v1, v2);
                city->addEdge(v2, v1);
            } else {
                city->addEdge(v1, v2);
            }
        }
        (city->irEVir()) ? cout << "1" : cout << "0";
        cout << endl;
        cin >> corners >> streets;
    }

    return 0;
}
