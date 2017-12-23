/*
 * This code is provied as a sample code of Hw 2 of "Theory of Computer Game".
 * The "genmove" function will randomly output one of the legal move.
 * This code can only be used within the class.
 *
 * 2015 Nov. Hung-Jui Chang
 * */
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>


#define BOARDSIZE        9
#define BOUNDARYSIZE    11
#define COMMANDLENGTH 1000
#define DEFAULTTIME     10
#define DEFAULTKOMI      7

#define MAXGAMELENGTH 1000
#define MAXSTRING       50
#define MAXDIRECTION     4

#define NUMINTERSECTION 81
#define HISTORYLENGTH   200

#define EMPTY            0
#define BLACK            1
#define WHITE            2
#define BOUNDARY         3

#define SELF             1
#define OPPONENT         2

#define NUMGTPCOMMANDS      15

#define LOCALVERSION      1
#define GTPVERSION        1
 
using namespace std;
int _board_size = BOARDSIZE;
int _board_boundary = BOUNDARYSIZE;
double _komi =  DEFAULTKOMI;
const int DirectionX[MAXDIRECTION] = {-1, 0, 1, 0};
const int DirectionY[MAXDIRECTION] = { 0, 1, 0,-1};
const char LabelX[]="0ABCDEFGHJ";

int change_turn[2] = {2,1};
char change_char_turn[2] = {'w','b'};

random_device rd;
mt19937 mt(rd());
uniform_int_distribution<int> dist(1, 100);

double final_score(int Board[BOUNDARYSIZE][BOUNDARYSIZE]);
void gtp_showboard(int Board[BOUNDARYSIZE][BOUNDARYSIZE]);


/*
 * This function reset the board, the board intersections are labeled with 0,
 * the boundary intersections are labeled with 3.
 * */
void reset(int Board[BOUNDARYSIZE][BOUNDARYSIZE]) {
    for (int i = 1 ; i <= BOARDSIZE; ++i) {
	for (int j = 1 ; j <= BOARDSIZE; ++j) {
	    Board[i][j] = EMPTY;
	}
    }
    for (int i = 0 ; i < BOUNDARYSIZE; ++i) {
	Board[0][i] = Board[BOUNDARYSIZE-1][i] = Board[i][0] = Board[i][BOUNDARYSIZE-1] = BOUNDARY;
    }
}

/*
 * This function return the total number of liberity of the string of (X, Y) and
 * the string will be label with 'label'.
 * */
int find_liberty(int X, int Y, int label, int Board[BOUNDARYSIZE][BOUNDARYSIZE], int ConnectBoard[BOUNDARYSIZE][BOUNDARYSIZE]) {
    // Label the current intersection
    ConnectBoard[X][Y] |= label;
    int total_liberty = 0;
    for (int d = 0 ; d < MAXDIRECTION; ++d) {
	// Check this intersection has been visited or not
	if ((ConnectBoard[X+DirectionX[d]][Y+DirectionY[d]] & (1<<label) )!= 0) continue;

	// Check this intersection is not visited yet
	ConnectBoard[X+DirectionX[d]][Y+DirectionY[d]] |=(1<<label);
	// This neighboorhood is empty
	if (Board[X+DirectionX[d]][Y+DirectionY[d]] == EMPTY){
	    total_liberty++;
	}
	// This neighboorhood is self stone
	else if (Board[X+DirectionX[d]][Y+DirectionY[d]] == Board[X][Y]) {
	    total_liberty += find_liberty(X+DirectionX[d], Y+DirectionY[d], label, Board, ConnectBoard);
	}
    }
    return total_liberty;
}

/*
 * This function count the liberties of the given intersection's neighboorhod
 * */
void count_liberty(int X, int Y, int Board[BOUNDARYSIZE][BOUNDARYSIZE], int Liberties[MAXDIRECTION]) {
    int ConnectBoard[BOUNDARYSIZE][BOUNDARYSIZE];
    // Initial the ConnectBoard
    for (int i = 0 ; i < BOUNDARYSIZE; ++i) {
	for (int j = 0 ; j < BOUNDARYSIZE; ++j) {
	    ConnectBoard[i][j] = 0;
	}
    }
    // Find the same connect component and its liberity
    for (int d = 0 ; d < MAXDIRECTION; ++d) {
	Liberties[d] = 0;
	if (Board[X+DirectionX[d]][Y+DirectionY[d]] == BLACK ||  
	    Board[X+DirectionX[d]][Y+DirectionY[d]] == WHITE    ) {
	    Liberties[d] = find_liberty(X+DirectionX[d], Y+DirectionY[d], d, Board, ConnectBoard);
	}
    }
}

/*
 * This function count the number of empty, self, opponent, and boundary intersections of the neighboorhod
 * and saves the type in NeighboorhoodState.
 * */
void count_neighboorhood_state(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn, int* empt, int* self, int* oppo ,int* boun, int NeighboorhoodState[MAXDIRECTION]) {
    for (int d = 0 ; d < MAXDIRECTION; ++d) {
	// check the number of nonempty neighbor
	switch(Board[X+DirectionX[d]][Y+DirectionY[d]]) {
	    case EMPTY:    (*empt)++; 
			   NeighboorhoodState[d] = EMPTY;
			   break;
	    case BLACK:    if (turn == BLACK) {
			       (*self)++;
			       NeighboorhoodState[d] = SELF;
			   }
			   else {
			       (*oppo)++;
			       NeighboorhoodState[d] = OPPONENT;
			   }
			   break;
	    case WHITE:    if (turn == WHITE) {
			       (*self)++;
			       NeighboorhoodState[d] = SELF;
			   }
			   else {
			       (*oppo)++;
			       NeighboorhoodState[d] = OPPONENT;
			   }
			   break;
	    case BOUNDARY: (*boun)++;
			   NeighboorhoodState[d] = BOUNDARY;
			   break;
	}
    }
}

/*
 * This function remove the connect component contains (X, Y) with color "turn" 
 * And return the number of remove stones.
 * */
int remove_piece(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn) {
    int remove_stones = (Board[X][Y]==EMPTY)?0:1;
    Board[X][Y] = EMPTY;
    for (int d = 0; d < MAXDIRECTION; ++d) {
	if (Board[X+DirectionX[d]][Y+DirectionY[d]] == turn) {
	    remove_stones += remove_piece(Board, X+DirectionX[d], Y+DirectionY[d], turn);
	}
    }
    return remove_stones;
}
/*
 * This function update Board with place turn's piece at (X,Y).
 * Note that this function will not check if (X, Y) is a legal move or not.
 * */
void update_board(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn) {
    int num_neighborhood_self = 0;
    int num_neighborhood_oppo = 0;
    int num_neighborhood_empt = 0;
    int num_neighborhood_boun = 0;
    int Liberties[4];
    int NeighboorhoodState[4];
    count_neighboorhood_state(Board, X, Y, turn,
	    &num_neighborhood_empt,
	    &num_neighborhood_self,
	    &num_neighborhood_oppo,
	    &num_neighborhood_boun, NeighboorhoodState);
    // check if there is opponent piece in the neighboorhood
    if (num_neighborhood_oppo != 0) {
	count_liberty(X, Y, Board, Liberties);
	for (int d = 0 ; d < MAXDIRECTION; ++d) {
	    // check if there is opponent component only one liberty
	    if (NeighboorhoodState[d] == OPPONENT && Liberties[d] == 1 && Board[X+DirectionX[d]][Y+DirectionY[d]]!=EMPTY) {
		remove_piece(Board, X+DirectionX[d], Y+DirectionY[d], Board[X+DirectionX[d]][Y+DirectionY[d]]);
	    }
	}
    }
    Board[X][Y] = turn;
}
/*
 * This function update Board with place turn's piece at (X,Y).
 * Note that this function will check if (X, Y) is a legal move or not.
 * */
int update_board_check(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn) {
    // Check the given coordination is legal or not
    if ( X < 1 || X > BOARDSIZE || Y < 1 || Y > BOARDSIZE || Board[X][Y]!=EMPTY)
	return 0;
    int num_neighborhood_self = 0;
    int num_neighborhood_oppo = 0;
    int num_neighborhood_empt = 0;
    int num_neighborhood_boun = 0;
    int Liberties[4];
    int NeighboorhoodState[4];
    count_neighboorhood_state(Board, X, Y, turn,
	    &num_neighborhood_empt,
	    &num_neighborhood_self,
	    &num_neighborhood_oppo,
	    &num_neighborhood_boun, NeighboorhoodState);
    // Check if the move is a legal move
    // Condition 1: there is a empty intersection in the neighboorhood
    int legal_flag = 0;
    count_liberty(X, Y, Board, Liberties);
    if (num_neighborhood_empt != 0) {
	legal_flag = 1;
    }
    else {
	// Condition 2: there is a self string has more than one liberty
	for (int d = 0; d < MAXDIRECTION; ++d) {
	    if (NeighboorhoodState[d] == SELF && Liberties[d] > 1) {
		legal_flag = 1;
	    }
	}
	if (legal_flag == 0) {
	// Condition 3: there is a opponent string has exactly one liberty
	    for (int d = 0; d < MAXDIRECTION; ++d) {
		if (NeighboorhoodState[d] == OPPONENT && Liberties[d] == 1) {
		    legal_flag = 1;
		}
	    }
	}
    }

    if (legal_flag == 1) {
    // check if there is opponent piece in the neighboorhood
	if (num_neighborhood_oppo != 0) {
	    for (int d = 0 ; d < MAXDIRECTION; ++d) {
		// check if there is opponent component only one liberty
		if (NeighboorhoodState[d] == OPPONENT && Liberties[d] == 1 && Board[X+DirectionX[d]][Y+DirectionY[d]]!=EMPTY) {
		    remove_piece(Board, X+DirectionX[d], Y+DirectionY[d], Board[X+DirectionX[d]][Y+DirectionY[d]]);
		}
	    }
	}
	Board[X][Y] = turn;
    }

    return (legal_flag==1)?1:0;
}

/*
 * This function return the number of legal moves with clor "turn" and
 * saves all legal moves in MoveList
 * */
int gen_legal_move(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int game_length, int GameRecord[BOUNDARYSIZE][BOUNDARYSIZE], int MoveList[HISTORYLENGTH]) {
    int NextBoard[BOUNDARYSIZE][BOUNDARYSIZE];
    int num_neighborhood_self = 0;
    int num_neighborhood_oppo = 0;
    int num_neighborhood_empt = 0;
    int num_neighborhood_boun = 0;
    int legal_moves = 0;
    int next_x, next_y;
    int Liberties[4];
    int NeighboorhoodState[4];
    bool eat_move = 0;
    for (int x = 1 ; x <= BOARDSIZE; ++x) 
	{
		for (int y = 1 ; y <= BOARDSIZE; ++y) 
		{
		    // check if current 
		    if (Board[x][y] == 0) 
			{
				// check the liberty of the neighborhood intersections
				num_neighborhood_self = 0;
				num_neighborhood_oppo = 0;
				num_neighborhood_empt = 0;
				num_neighborhood_boun = 0;
				// count the number of empy, self, opponent, and boundary neighboorhood
				count_neighboorhood_state(Board, x, y, turn,
					&num_neighborhood_empt,
					&num_neighborhood_self,
					&num_neighborhood_oppo,
					&num_neighborhood_boun, NeighboorhoodState);
				// check if the emtpy intersection is a legal move
				next_x = next_y = 0;
				eat_move = 0;
				count_liberty(x, y, Board, Liberties);
				int recopy = 0;
				// Case 1: exist empty intersection in the neighborhood
				 if (num_neighborhood_empt > 0) 
				 {
				     next_x = x;
				     next_y = y;
				     // check if it is a capture move
				     for (int d = 0 ; d < MAXDIRECTION; ++d) 
					 {
						 if (NeighboorhoodState[d] == OPPONENT && Liberties[d] == 1) 
						 {
						 	if( GameRecord[x+DirectionX[d]][y+DirectionY[d]] != game_length )
						 	{
						 		//cerr << "Test1" << endl;
						 		eat_move = 1;
							}
							else
							{
								int test = 0;
								for (int e = 0 ; e < MAXDIRECTION; ++e) 
								{
									int move_x = x+DirectionX[d];
									int move_y = y+DirectionY[d];
									if( move_x+DirectionX[e] == -1 || move_x+DirectionX[e] == 11 || move_y+DirectionY[e]== -1 || move_y+DirectionY[e] == 11 )
										continue;
									else
									{
										if( Board[move_x][move_y] == Board[move_x+DirectionX[e]][move_y+DirectionY[e]] )
										{
											test = 1;
											break;
										}
									}
								}
								if(test==0)
								{
									recopy = 1;
								}
								else
								{
									eat_move = 1;
								}
							}
						 }
				     }
		
				 }
				 // Case 2: no empty intersection in the neighborhood
				 else {
				    // Case 2.1: Surround by the self piece
				    if (num_neighborhood_self + num_neighborhood_boun == MAXDIRECTION) {
					int check_flag = 0, check_eye_flag = num_neighborhood_boun;
					for (int d = 0 ; d < MAXDIRECTION; ++d) {
					    // Avoid fill self eye
					    if (NeighboorhoodState[d]==SELF && Liberties[d] > 1) {
						check_eye_flag++;
					    }
					    // Check if there is one self component which has more than one liberty
					    if (NeighboorhoodState[d]==SELF && Liberties[d] > 1) {
						check_flag = 1;
					    }
					}
					if (check_flag == 1 && check_eye_flag!=4) {
					    next_x = x;
					    next_y = y;
					}
				    }	
				    // Case 2.2: Surround by opponent or both side's pieces.
				    else if (num_neighborhood_oppo > 0) {
					int check_flag = 0;
					int eat_flag = 0;
					int temp = 0;
					for (int d = 0 ; d < MAXDIRECTION; ++d) {
					    // Check if there is one self component which has more than one liberty
					    if (NeighboorhoodState[d]==SELF && Liberties[d] > 1) {
						check_flag = 1;
					    }
					    // Check if there is one opponent's component which has exact one liberty
					    if (NeighboorhoodState[d]==OPPONENT && Liberties[d] == 1) {
					    	temp = d;
							eat_flag = 1;
					    }
					}
					if (check_flag == 1) {
					    next_x = x;
					    next_y = y;
					    if (eat_flag == 1) {
						 		eat_move = 1;
					    }
					}
					else { // check_flag == 0
					    if (eat_flag == 1) 
						{
							next_x = x;
							next_y = y;
							if( GameRecord[x+DirectionX[temp]][y+DirectionY[temp]] != game_length )
						 	{
							 	//cout << x+DirectionX[temp] << y+DirectionY[temp] << " " << game_length << endl;
							 	//recopy = 0;
							 	eat_move = 1;
							}
							else
							{
								int test = 0;
								for (int d = 0 ; d < MAXDIRECTION; ++d) 
								{
									int move_x = x+DirectionX[temp];
									int move_y = y+DirectionY[temp];
									if( move_x+DirectionX[d] == -1 || move_x+DirectionX[d] == 11 || move_y+DirectionY[d]== -1 || move_y+DirectionY[d] == 11 )
										continue;
									else
									{
										if( Board[move_x][move_y] == Board[move_x+DirectionX[d]][move_y+DirectionY[d]] )
										{
											test = 1;
											break;
										}
									}
								}
								if(test==0)
								{
									recopy = 1;
								}
								else
								{
									eat_move = 1;
								}
							}

					    }
					}
				    }	
				 }
				 if (next_x !=0 && next_y !=0 && recopy!=1) {
				 // copy the current board to next board
				    for (int i = 0 ; i < BOUNDARYSIZE; ++i) 
					{
						for (int j = 0 ; j < BOUNDARYSIZE; ++j) 
						{
						    NextBoard[i][j] = Board[i][j];
						}
				    }
				    // do the move
				    // The move is a capture move and the board needs to be updated.
				    if (eat_move == 1) {
						update_board(NextBoard, next_x, next_y, turn);
				    }
				    else {
					NextBoard[x][y] = turn;
				    }
				    MoveList[legal_moves] = eat_move * 100 + next_x * 10 + y ;
					legal_moves++;
				 }
		    }
		}
    }
    return legal_moves;
}

/*
 * This function randomly selects one move from the MoveList.
 * */
int rand_pick_move(int num_legal_moves, int MoveList[HISTORYLENGTH]) {
    if (num_legal_moves == 0)
	return 0;
    else {
	int move_id = dist(mt)%num_legal_moves;
	return MoveList[move_id];
    }
}
/*
 * This function update the Board with put 'turn' at (x,y)
 * where x = (move % 100) / 10 and y = move % 10.
 * Note this function will not check 'move' is legal or not.
 * */
void do_move(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int move) {
    int move_x = (move % 100) / 10;
    int move_y = move % 10;
    if (move<100) {
	Board[move_x][move_y] = turn;
    }
    else {
	update_board(Board, move_x, move_y, turn);
    }

}

/* 
 * This function records the current game baord with current
 * game length "game_length"
 * */
void record(int GameRecord[BOUNDARYSIZE][BOUNDARYSIZE], int move, int game_length) 
{
	int move_i = (move%100)/10;
	int move_j = (move%10);
	GameRecord[move_i][move_j] = game_length;
}

struct tree                       
{
   int data[BOUNDARYSIZE][BOUNDARYSIZE][2];
   int Board[BOUNDARYSIZE][BOUNDARYSIZE];
   int GameRecord[BOUNDARYSIZE][BOUNDARYSIZE];
   int count;
   int path;
   int depth;
   int best_move;
   vector <tree* > son;
   tree *parent;
   
};

tree* GetNewNode(tree* root, int Board[BOUNDARYSIZE][BOUNDARYSIZE], int GameRecord[BOUNDARYSIZE][BOUNDARYSIZE], int game_length, int path)
{
	tree* newNode = new tree();
	memcpy(newNode->Board, Board, sizeof(newNode->Board));
	memcpy(newNode->GameRecord, GameRecord, sizeof(newNode->GameRecord));
	newNode->depth = game_length;
	newNode->path = path;
	newNode->count = 0;
	newNode->best_move = 0;
	newNode->parent = root;
	return newNode;
}

tree* UCT_Selection(tree* root, int turn, int game_length)
{
	//cout << "UCT_Selection" << endl;
	if(root == NULL)
	{
		return root;
	}
	double max = 0;

	double c = sqrt(2);
	int best_move = 0;
	int MoveList[HISTORYLENGTH];
	
	int num_legal_moves = gen_legal_move(root->Board, turn, game_length, root->GameRecord, MoveList);
	//gtp_showboard(root->Board);
	for(int i=0; i<num_legal_moves; i++)
	{
		int move_x = (MoveList[i] % 100) / 10;
	    int move_y = MoveList[i] % 10;
	    double win_rate;
		double exploration;
		if(root->data[move_x][move_y][root->depth%2] == -1)
		{
			continue;
		}
		if(root->data[move_x][move_y][0]+root->data[move_x][move_y][1] == 0)
		{
			continue;
		}
		else
		{
			win_rate = (double) root->data[move_x][move_y][root->depth%2] / (double)(root->data[move_x][move_y][0]+root->data[move_x][move_y][1]);
			exploration = sqrt( (double) log(root->count) / (double)(root->data[move_x][move_y][0]+root->data[move_x][move_y][1]) );
			//exploration = 0;
		}
		if( max <  win_rate+c*exploration && root->Board[move_x][move_y] == EMPTY )
		{
			max = win_rate+ c*exploration;
			best_move = move_x*10+move_y;
			//cout << move_x << move_x << ",max=" << max << endl;
		}
	}
	
	root->best_move = best_move;
	if(root->son.size() != 0)
	{
		//cout << "秈ㄓ" << endl;
		tree* temp;
		for(int i=0;i<root->son.size();i++)
		{
			if( root->son[i]->path == root->best_move )
			{
				return UCT_Selection(root->son[i], change_turn[(root->depth+turn)%2], game_length+1);
			}
		}
	}
	
	return root;
}

tree* UCT_Expansion(tree* root, int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int game_length, int best_move)
{
	//cout << "UCT_Expansion" << endl;
	if(root != NULL)
	{
		//gtp_showboard(root->Board);
		//cout << "root->son.size()= " << root->son.size() << endl;
		root->son.push_back(GetNewNode(root,root->Board,root->GameRecord,root->depth+1,best_move));
		update_board(root->son[root->son.size()-1]->Board, (best_move%100)/10, best_move%10, change_turn[(game_length+turn)%2]);
		//gtp_showboard(root->son[root->son.size()-1]->Board);
		return root->son[root->son.size()-1];
	}
}

void UCT_Backward_propagation(tree* root)
{
	//cout << "UCT_Backward_propagation" << endl;
	if( root == NULL )
	{
	}
	else
	{
		if( root->parent != NULL)
		{
			tree* root_parent;
			tree* root_now = root;
			while( root_now->parent != NULL)
			{
				//cout << "root_now->path=" << root_now->path << endl;
				root_parent = root_now->parent;
				int move_x = (root_now->path % 100) / 10;
				int move_y = root_now->path % 10;   
				for(int i=1; i<BOARDSIZE; i++)
				{
					for(int j=1; j<BOARDSIZE; j++)
					{
						root_parent->data[move_x][move_y][0] += root->data[i][j][0];
						root_parent->data[move_x][move_y][1] += root->data[i][j][1];
					}
				}
				root_parent->data[0][0][0] += root->data[0][0][0];
				root_parent->data[0][0][1] += root->data[0][0][1];
				root_parent->count += root->count;
				root_now = root_parent;
			}
		}
		
	}
}

void Progressive_pruning(tree* root)
{
	//cout << "Progressive_pruning" << endl;
	double sum = 0;
	double max = 0;
	double count = 0;
	for(int i=1;i<=BOARDSIZE;i++)
	{
		for(int j=1;j<=BOARDSIZE;j++)
		{
			//cerr << root->data[i][j][root->depth%2] << " ";
			if( max < root->data[i][j][root->depth%2] )
			{
				max = root->data[i][j][root->depth%2];
			}
			if(root->data[i][j][0]+root->data[i][j][1] != 0)
			{
				sum += (double) root->data[i][j][root->depth%2];
				count ++;
			}
		}
		//cerr << endl;
	}
	sum /= count;
	//cerr << "sum=" << sum << endl;
	double SD = 0;
	for(int i=1;i<=BOARDSIZE;i++)
	{
		for(int j=1;j<=BOARDSIZE;j++)
		{
			if(root->data[i][j][0]+root->data[i][j][1] != 0)
				SD += (sum-(double)root->data[i][j][root->depth%2])*(sum-(double)root->data[i][j][root->depth%2]);
		}
	}
	SD /= count;
	SD = sqrt(SD);
	//cerr << "SD=" << SD << endl;
	for(int i=1;i<=BOARDSIZE;i++)
	{
		for(int j=1;j<=BOARDSIZE;j++)
		{
			if( max - SD > root->data[i][j][root->depth%2] + SD )
			{
				//cout << i << j << "砆奔" << endl;
				root->data[i][j][root->depth%2]= -1;
			}
		}
	}
	
}

void Print(tree* root)
{
	if(root != NULL)
	{
		if(root->son.size()!=0)
		{
			cerr << "瞏" << root->depth << " " ;
			cerr << "parent=" << root->path << " "; 
			for(int i=0;i<root->son.size();i++)
			{
				cerr << root->son[i]->path << " ";
			}
			cerr << endl;
			for(int i=0;i<root->son.size();i++)
			{
				Print(root->son[i]);
			}
		}
	}
}

void Print_PV(tree* root,int turn)
{
	if(root != NULL)
	{
		//gtp_showboard(root->Board);
		double c = sqrt(2);
		double max = 0;
		double win_rate = 0;
		double exploration = 0;
		int best_move = 0;
		for(int i=1; i<=BOARDSIZE; i++)
		{
			for(int j=1; j<=BOARDSIZE; j++)
			{
				if(root->data[i][j][root->depth%2] == -1 || root->data[i][j][root->depth%2] == 0)
				{
					continue;
				}
				win_rate = (double) root->data[i][j][root->depth%2] / (double)(root->data[i][j][0]+root->data[i][j][1]);
				//exploration = sqrt( (double) log(root->count) / (double)(root->data[i][j][0]+root->data[i][j][1]) );
				exploration = 0;
				if( max <  win_rate+c*exploration && root->Board[i][j] == EMPTY )
				{
					max = win_rate+c*exploration;
					best_move = i*10+j;
				}
			}
		}
		int move_x = (best_move % 100) / 10;
		int move_y = best_move % 10;
		cerr << root->depth+1 << ". " << change_char_turn[turn%2] << " " << LabelX[move_y] << 10-move_x;
		max = (double)root->data[move_x][move_y][root->depth%2]/(double)(root->data[move_x][move_y][0]+root->data[move_x][move_y][1]);
		if(root->depth%2==0)
			cerr << ", win rate=" << max << ", sim=" << root->data[move_x][move_y][0]+root->data[move_x][move_y][1] << endl;
		else
			cerr << ", win rate=" << (double)1-max << ", sim=" << root->data[move_x][move_y][0]+root->data[move_x][move_y][1] << endl;
		
		//cout << "best_move=" << best_move << endl;
		for(int i=0; i<root->son.size(); i++)
		{
			//cout << "son[" << i << "]=" << root->son[i]->path << endl;
			if(root->son[i]->path == best_move )
			{
				if(turn == 1)
					Print_PV(root->son[i],2);
				else
					Print_PV(root->son[i],1);
			}
		}
	}
}

int OffLine(int Board[BOUNDARYSIZE][BOUNDARYSIZE],int game_length)
{
	if(game_length == 0)
	{
		return 55;
	}
	else if(game_length == 1)
	{
		if(Board[3][5] == EMPTY)
			return 35;
		else if(Board[5][3] == EMPTY)
			return 53;
	}
	else if(game_length == 2)
	{
		if(Board[3][5] == EMPTY)
			return 35;
		else if(Board[5][3] == EMPTY)
			return 53;
	}
	else
	{
		if(Board[3][5] == EMPTY)
			return 35;
		else if(Board[5][3] == EMPTY)
			return 53;
		else if(Board[5][7] == EMPTY)
			return 57;
		else
			return 75;
	}
}

/* 
 * This function randomly generate one legal move (x, y) with return value x*10+y,
 * if there is no legal move the function will return 0.
 * */
int genmove(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int time_limit, int game_length, int GameRecord[BOUNDARYSIZE][BOUNDARYSIZE]) {

	
	if(game_length < 4)
    {
    	int move = OffLine(Board,game_length);
    	do_move(Board, turn, move);
    	return move;
	}
	
	clock_t start_t, end_t, now_t;
    // record start time
    start_t = clock();
    // calculate the time bound
    end_t = start_t + CLOCKS_PER_SEC * time_limit;

    int MoveList[HISTORYLENGTH];
    int num_legal_moves = 0;
    int return_move = 0;

	int Simulation_Board[BOUNDARYSIZE][BOUNDARYSIZE];
	int Simulation_Board_GameRecord[BOUNDARYSIZE][BOUNDARYSIZE];
	
	//cout << "game_length=" << game_length << endl;
	num_legal_moves = gen_legal_move(Board, turn, game_length , GameRecord, MoveList);
    if(num_legal_moves == 0)
	{
		//cerr << "礚" << endl;
		return 0;
	}
	else
	{
		/*for(int i=0; i<num_legal_moves; i++)
	    {
	    	cerr << MoveList[i] << " ";
		}
		cerr << endl;*/
	}
	/*cerr << endl;
	for(int i=1; i<=BOARDSIZE; i++)
	{
		for(int j=1; j<=BOARDSIZE; j++)
		{
			cerr << GameRecord[i][j] << " ";
		}
		cerr << endl;
	}
	cerr << endl;
	return 0;*/
	
	tree* root;
	tree* PV;
	PV = GetNewNode(root,Board,GameRecord,0,0);
	root = PV;
	int count = 0;
	
	while(true)
	{
		if( end_t < clock() )
		{
			break;
		}
		/*if(count>1)
		{
			break;
		}*/

		if(count != 0)
		{
			PV = UCT_Selection(root,turn,game_length);
			num_legal_moves = gen_legal_move(PV->Board, change_turn[(PV->depth+turn)%2], game_length+PV->depth+1 , PV->GameRecord, MoveList);
			//gtp_showboard(PV->Board);
			//cerr << "count=" << count << ",猭˙计:" << num_legal_moves << endl;
			if(num_legal_moves != 0)
			{
				PV = UCT_Expansion(PV,PV->Board,turn,PV->depth,PV->best_move);
				record(PV->GameRecord,PV->path,game_length+PV->depth);
			}
			else
			{
				record(PV->GameRecord,0,game_length+PV->depth);
				num_legal_moves = gen_legal_move(PV->Board, change_turn[(PV->depth+turn+1)%2], game_length+PV->depth+2 , PV->GameRecord, MoveList);
				if(num_legal_moves == 0)
				{
					double result;
					int win;
					result = final_score(PV->Board);
					//cerr << "AAA**********FINAL_BOARD**********" << endl;
					//gtp_showboard(PV->Board);
					result -= _komi;
					if (result >= 0.0) 
					{ // Black win
						win = 1;
					}
					if (result < 0.0) 
					{ // White win
						win = 2;
					}
					if( win == turn )
					    	PV->data[(PV->path%100)/10][PV->path%10][0]++;
					    else
					    	PV->data[(PV->path%100)/10][PV->path%10][1]++;
					PV->count++;
					//cerr << "win=" << win << endl;
					//cerr << "turn=" << change_turn[(PV->depth+turn)%2] << endl;
				}
			}
		}
		//cerr << "************ROOT STRUCTURE***************" << endl;
		//Print(root);
		//cerr << endl;
		
		//cout << "depth= " << PV->depth << endl;
		num_legal_moves = gen_legal_move(PV->Board, change_turn[(PV->depth+turn)%2], game_length+PV->depth+1 , PV->GameRecord, MoveList);
    	//cerr << "猭˙计:" << num_legal_moves << endl;
    	


		for(int i=0; i < num_legal_moves ; i++)
	    {
	    	//cerr << "i=" << i << endl;
	    	for(int j=0; j<10; j++)
	    	{
	    		if( end_t < clock() )
				{
					break;
				}
				PV->count += 10;
				
		    	int move_x = (MoveList[i] % 100) / 10;
			    int move_y = MoveList[i] % 10;
			    
			    /*if(PV->data[move_x][move_y][PV->depth%2] == -1)
			    {
			    	break;
				}*/
			    
			    if( (PV->depth+game_length) < 14)
				{
					if( move_x < 3 || move_y < 3 || move_x > 7 || move_y > 7)
						continue;
				}
				if( (PV->depth+game_length) < 40)
				{
					if( move_x < 2 || move_y < 2 || move_x > 8 || move_y > 8)
						continue;
				}
				
				//cerr << "move:" << MoveList[i] << endl;
			    memcpy(Simulation_Board, PV->Board, sizeof(Simulation_Board));
			    memcpy(Simulation_Board_GameRecord, GameRecord, sizeof(Simulation_Board_GameRecord));
			    if ( MoveList[i] < 100 )
				{
					Simulation_Board[move_x][move_y] = change_turn[(game_length+PV->depth+1)%2];
			    }
			    else 
				{
					update_board(Simulation_Board, move_x, move_y, change_turn[(game_length+PV->depth+1)%2]);
			    }
			    record(Simulation_Board_GameRecord,MoveList[i],game_length+PV->depth+1);
			    //gtp_showboard(Simulation_Board);

			    
			    int win = 0;
			    int length = game_length+PV->depth+1;
			    int temp = length;
			    int pass = 0;
			    while(true)
			    {
			    	if( end_t < clock() )
					{
						break;
					}
					/*if(length < temp+12)
					{
						gtp_showboard(Simulation_Board);
						cerr << endl;
						for(int i=1; i<=BOARDSIZE; i++)
						{
							for(int j=1; j<=BOARDSIZE; j++)
							{
								cerr << Simulation_Board_GameRecord[i][j] << " ";
							}
							cerr << endl;
						}
						cerr << endl;
					}*/
			    	//cerr << "length = " << length << endl;
			    	//gtp_showboard(Simulation_Board);
			    	int num_legal_moves = 0;
					int MoveList[HISTORYLENGTH];
			    	num_legal_moves = gen_legal_move(Simulation_Board, change_turn[(length+1)%2] , length, Simulation_Board_GameRecord, MoveList);
					
					/*cerr << change_turn[(length+1)%2] << "家览猭˙计:" << num_legal_moves << endl;
					cerr << change_turn[(length+1)%2] << ": " ;
					for(int i=0; i<num_legal_moves; i++)
				    {
				    	cerr << MoveList[i] << " ";
					}
					cerr << endl;*/
					
					if(num_legal_moves == 0)
					{
						pass++;
					}
					else
					{
						pass = 0;
					}
					
					if(pass == 2)
					{
						double result;
				    	result = final_score(Simulation_Board);
				    	result -= _komi;
					    if (result >= 0.0) 
						{ // Black win
							//cout << "B+" << result << endl << endl<< endl;
							win = 1;
							break;
					    }
					    if (result < 0.0) 
						{ // White win
							//cout << "W+" << -result << endl << endl<< endl;
							win = 2;
							break;
					    }
					}
					else
					{
						int return_move = rand_pick_move(num_legal_moves, MoveList);
						//cerr << "move " << change_turn[(turn+length)%2] << ": " << return_move << endl;

						int move_x = (return_move % 100) / 10;
				    	int move_y = return_move % 10;
						if ( return_move < 100 ) 
						{
							Simulation_Board[move_x][move_y] = change_turn[(length+1)%2];
					    }
					    else 
						{
							update_board(Simulation_Board, move_x, move_y, change_turn[(length+1)%2]);
					    }
					    record(Simulation_Board_GameRecord,return_move,length+1);
						//Simulation_Board[move_x][move_y] = turn;
						//update_board(Simulation_Board, move_x, move_y, turn);
						//gtp_showboard(Simulation_Board);
						
						length++;
					}
				}
			    if( win == turn )
			    	PV->data[move_x][move_y][0]++;
			    else
			    	PV->data[move_x][move_y][1]++;
			}
		    
		}
		//cout << "PV->count=" << PV->count << endl;
		count++;
		UCT_Backward_propagation(PV);
		Progressive_pruning(PV);
	}
	//cout << (double)(clock()-start_t)/(double)CLOCKS_PER_SEC << endl;
	
	/*cerr << endl;
	for(int i=0; i<BOUNDARYSIZE; i++)
	{
		for(int j=0; j<BOUNDARYSIZE; j++)
		{
			cerr << (double)root->data[i][j][0]/(double)(root->data[i][j][0]+root->data[i][j][1]) << "\t";
		}
		cerr << endl;
	}
	cerr << endl;
	

	//gtp_showboard(root->son[2]->Board);
	cerr << endl;
	for(int i=0; i<BOUNDARYSIZE; i++)
	{
		for(int j=0; j<BOUNDARYSIZE; j++)
		{
			cerr << root->data[i][j][0] << " ";
		}
		cerr << endl;
	}
	cerr << endl;
	for(int i=0; i<BOUNDARYSIZE; i++)
	{
		for(int j=0; j<BOUNDARYSIZE; j++)
		{
			cerr << root->data[i][j][1] << " ";
		}
		cerr << endl;
	}
	cerr << endl;*/
	
	//return 0;
	
	num_legal_moves = gen_legal_move(Board, turn, game_length, GameRecord, MoveList);
	double max = 0;
	int best_move = 0;
	
	for(int i=0; i<num_legal_moves; i++)
	{
		int move_x = (MoveList[i] % 100) / 10;
		int move_y = MoveList[i] % 10;
		if( root->data[move_x][move_y][0] == -1 )
		{
			continue;
		}
		double win_rate = (double)root->data[move_x][move_y][0]/(double)(root->data[move_x][move_y][0]+root->data[move_x][move_y][1]);
		if( max < win_rate && Board[move_x][move_y] == EMPTY )
		{
			max = win_rate;
			best_move = 10*move_x+move_y;
		}
	}
	if(max < 0.15 && num_legal_moves < 4 )
	{
		return 0;
	}
	//cout << "best_move= " << best_move << endl;
	//cout << "win_rate= " << max << endl;

	//cout << (double)(clock()-start_t)/(double)CLOCKS_PER_SEC << endl;
	
	Print_PV(root,turn);
	
    //do_move(Board, turn, best_move);
    update_board(Board, (best_move%100)/10, best_move%10, turn);
    return best_move % 100;
}

/*
 * This function counts the number of points remains in the board by Black's view
 * */
double final_score(int Board[BOUNDARYSIZE][BOUNDARYSIZE]) {
    int black, white;
    black = white = 0;
    int is_black, is_white;
    for (int i = 1 ; i <= BOARDSIZE; ++i) {
	for (int j = 1; j <= BOARDSIZE; ++j) {
	    switch(Board[i][j]) {
		case EMPTY:
		    is_black = is_white = 0;
		    for(int d = 0 ; d < MAXDIRECTION; ++d) {
			if (Board[i+DirectionX[d]][j+DirectionY[d]] == BLACK) is_black = 1;
			if (Board[i+DirectionX[d]][j+DirectionY[d]] == WHITE) is_white = 1;
		    }
		    if (is_black + is_white == 1) {
			black += is_black;
			white += is_white;
		    }
		    break;
		case WHITE:
		    white++;
		    break;
		case BLACK:
		    black++;
		    break;
	    }
	}
    }
    return black - white;
}
/* 
 * Following are commands for Go Text Protocol (GTP)
 *
 * */
const char *KnownCommands[]={
    "protocol_version",
    "name",
    "version",
    "known_command",
    "list_commands",
    "quit",
    "boardsize",
    "clear_board",
    "komi",
    "play",
    "genmove",
    "undo",
    "quit",
    "showboard",
    "final_score"
};

void gtp_final_score(int Board[BOUNDARYSIZE][BOUNDARYSIZE]) {
    double result;
    result = final_score(Board);
    result -= _komi;
    cout << "= ";
    if (result > 0.0) { // Black win
	cout << "B+" << result << endl << endl<< endl;;
    }
    if (result < 0.0) { // White win
	cout << "W+" << -result << endl << endl<< endl;;
    }
    else { // draw
	cout << "0" << endl << endl<< endl;;
    }
}
void gtp_undo(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int game_length, int GameRecord[BOUNDARYSIZE][BOUNDARYSIZE]) {
    if (game_length!=0) 
	{
		for (int i = 1; i <= BOARDSIZE; ++i) 
		{
		    for (int j = 1; j <= BOARDSIZE; ++j) 
			{
				Board[i][j] = GameRecord[i][j];
		    }
		}
    }
    cout << "= " << endl << endl;
}
void gtp_showboard(int Board[BOUNDARYSIZE][BOUNDARYSIZE]) {
    for (int i = 1; i <=BOARDSIZE; ++i) {
	cerr << "#";
	cerr <<10-i;
	for (int j = 1; j <=BOARDSIZE; ++j) {
	    switch(Board[i][j]) {
		case EMPTY: cerr << " .";break;
		case BLACK: cerr << " X";break;
		case WHITE: cerr << " O";break;
	    }
	}
	cerr << endl;
    }
    cerr << "#  ";
    for (int i = 1; i <=BOARDSIZE; ++i) 
	cerr << LabelX[i] <<" ";
    cerr << endl;
    cerr << endl;

}
void gtp_protocol_version() {
    cout <<"= 2"<<endl<< endl;
}
void gtp_name() {
    cout <<"= TCG-randomGo99" << endl<< endl;
}
void gtp_version() {
    cout << "= 1.02" << endl << endl;
}
void gtp_list_commands(){
    cout <<"= ";
    for (int i = 0 ; i < NUMGTPCOMMANDS; ++i) {
	cout <<KnownCommands[i] << endl;
    }
    cout << endl;
}
void gtp_known_command(const char Input[]) {
    for (int i = 0 ; i < NUMGTPCOMMANDS; ++i) {
	if (strcmp(Input, KnownCommands[i])==0) {
	    cout << "= true" << endl<< endl;
	    return;
	}
    }
    cout << "= false" << endl<< endl;
}
void gtp_boardsize(int size) {
    if (size!=9) {
	cout << "? unacceptable size" << endl<< endl;
    }
    else {
	_board_size = size;
	cout << "= "<<endl<<endl;
    }
}
void gtp_clear_board(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int NumCapture[]) {
    reset(Board);
    NumCapture[BLACK] = NumCapture[WHITE] = 0;
    cout << "= "<<endl<<endl;
}
void gtp_komi(double komi) {
    _komi = komi;
    cout << "= "<<endl<<endl;
}
void gtp_play(char Color[], char Move[], int Board[BOUNDARYSIZE][BOUNDARYSIZE], int game_length, int GameRecord[BOUNDARYSIZE][BOUNDARYSIZE]) {
    int turn, move_i, move_j;
    if (Color[0] =='b' || Color[0] == 'B')
	turn = BLACK;
    else
	turn = WHITE;
    if (strcmp(Move, "PASS") == 0 || strcmp(Move, "pass")==0) 
	{
    }
    else {
	// [ABCDEFGHJ][1-9], there is no I in the index.
	Move[0] = toupper(Move[0]);
	move_j = Move[0]-'A'+1;
	if (move_j == 10) move_j = 9;
	move_i = 10-(Move[1]-'0');
	update_board(Board, move_i, move_j, turn);
	record(GameRecord, 10*move_i+move_j, game_length+1);
    }
    cout << "= "<<endl<<endl;
}
void gtp_genmove(int Board[BOUNDARYSIZE][BOUNDARYSIZE], char Color[], int time_limit, int game_length, int GameRecord[BOUNDARYSIZE][BOUNDARYSIZE]){
    int turn = (Color[0]=='b'||Color[0]=='B')?BLACK:WHITE;
    int move = genmove(Board, turn, time_limit, game_length, GameRecord);
    int move_i, move_j;
    record(GameRecord, move, game_length+1);
    /*cerr << endl;
	for(int i=0; i<BOUNDARYSIZE; i++)
	{
		for(int j=0; j<BOUNDARYSIZE; j++)
		{
			cerr << GameRecord[i][j] << " ";
		}
		cerr << endl;
	}
	cerr << endl;*/
    if (move==0) 
	{
		cout << "= PASS" << endl<< endl<< endl;
    }
    else 
	{
		move_i = (move%100)/10;
		move_j = (move%10);
		cerr << "#turn("<<game_length<<"): (move, move_i,move_j)" << turn << ": " << move<< " " << move_i << " " << move_j << endl;
		cout << "= " << LabelX[move_j]<<10-move_i<<endl<< endl;
    }
}
/*
 * This main function is used of the gtp protocol
 * */
void gtp_main(int display) {
    char Input[COMMANDLENGTH]="";
    char Command[COMMANDLENGTH]="";
    char Parameter[COMMANDLENGTH]="";
    char Move[4]="";
    char Color[6]="";
    int ivalue;
    double dvalue;
    int Board[BOUNDARYSIZE][BOUNDARYSIZE]={{0}};
    int NumCapture[3]={0};// 1:Black, 2: White
    int time_limit = DEFAULTTIME;
    int GameRecord[BOUNDARYSIZE][BOUNDARYSIZE]={{0}};
    int game_length = 0;
    if (display==1) {
	gtp_list_commands();
	gtp_showboard(Board);
    }
    while (gets(Input) != 0) {
	sscanf(Input, "%s", Command);
	if (Command[0]== '#')
	    continue;

	if (strcmp(Command, "protocol_version")==0) {
	    gtp_protocol_version();
	}
	else if (strcmp(Command, "name")==0) {
	    gtp_name();
	}
	else if (strcmp(Command, "version")==0) {
	    gtp_version();
	}
	else if (strcmp(Command, "list_commands")==0) {
	    gtp_list_commands();
	}
	else if (strcmp(Command, "known_command")==0) {
	    sscanf(Input, "known_command %s", Parameter);
	    gtp_known_command(Parameter);
	}
	else if (strcmp(Command, "boardsize")==0) {
	    sscanf(Input, "boardsize %d", &ivalue);
	    gtp_boardsize(ivalue);
	}
	else if (strcmp(Command, "clear_board")==0) {
	    gtp_clear_board(Board, NumCapture);
	    game_length = 0;
	}
	else if (strcmp(Command, "komi")==0) {
	    sscanf(Input, "komi %lf", &dvalue);
	    gtp_komi(dvalue);
	}
	else if (strcmp(Command, "play")==0) {
	    sscanf(Input, "play %s %s", Color, Move);
	    gtp_play(Color, Move, Board, game_length, GameRecord);
	    game_length++;
	    if (display==1) {
		gtp_showboard(Board);
	    }
	}
	else if (strcmp(Command, "genmove")==0) {
	    sscanf(Input, "genmove %s", Color);
	    gtp_genmove(Board, Color, time_limit, game_length, GameRecord);
	    game_length++;
	    if (display==1) {
		gtp_showboard(Board);
	    }
	}
	else if (strcmp(Command, "quit")==0) {
	    break;
	}
	else if (strcmp(Command, "showboard")==0) {
	    gtp_showboard(Board);
	}
	else if (strcmp(Command, "undo")==0) {
	    game_length--;
	    gtp_undo(Board, game_length, GameRecord);
	    if (display==1) {
		gtp_showboard(Board);
	    }
	}
	else if (strcmp(Command, "final_score")==0) {
	    if (display==1) {
		gtp_showboard(Board);
	    }
	    gtp_final_score(Board);
	}
    }
}
int main(int argc, char* argv[]) {
//    int type = GTPVERSION;// 1: local version, 2: gtp version
    int type = GTPVERSION;// 1: local version, 2: gtp version
    int display = 0; // 1: display, 2 nodisplay
    if (argc > 1) {
	if (strcmp(argv[1], "-display")==0) {
	    display = 1;
	}
	if (strcmp(argv[1], "-nodisplay")==0) {
	    display = 0;
	}
    }
    gtp_main(display);
    return 0;
}
